"""Parse output of TERMinator into :code:`.etab` files for use in MST.

Usage:
    .. code-block::

        python to_etab.py \\
            --output_dir <folder_with_net.out> \\
            --dtermen_data <dtermen_data_root> \\
            --num_cores <num_processes> \\
            [-u]

See :code:`python to_etab.py --help` for more info.
"""
import argparse
import multiprocessing as mp
import os
import pickle
from selectors import EpollSelector
import sys
import time
import traceback
import glob
import re
import numpy as np
from tqdm import tqdm

from terminator.utils.common import int_to_3lt_AA,AA_to_int,aa_to_similar

sys.path.append("../preprocessing")
from scripts.data.preprocessing.parseEtab import parseEtab
# pylint: disable=wrong-import-position,wrong-import-order,redefined-outer-name,unspecified-encoding

# for autosummary import purposes
sys.path.insert(0, os.path.dirname(__file__))
from search_utils import find_dtermen_folder


# print to stderr
def eprint(*args, **kwargs):
    """Print to stderr rather than stdout"""
    print(*args, file=sys.stderr, **kwargs)


# pylint: disable=broad-except
def _to_etab_file_wrapper(etab_matrix, E_idx, idx_dict, out_path):
    """Wrapper for _to_etab_file that does error handling"""
    try:
        return to_etab_file(etab_matrix, E_idx, idx_dict, out_path)
    except Exception:
        eprint(out_path)
        eprint(idx_dict)
        traceback.print_exc()
        return False, out_path

# Creates a new etab from an old one, where only positions matching selChainID are included and the rest are assumed to be the sequence in the PDB file.
def select_variable_positions(etab_matrix, E_idx, idx_dict, seq_dict, sel_chain_id):
    """Creates a new etab from an old one, where only positions matching selChainID are included and the rest are assumed to be the sequence in the PDB file.

    Args
    ====
    etab_matrix : np.ndarray
        Etab outputted by TERMinator
    E_idx : np.ndarray
        Indexing matrix associated with :code:`etab_matrix`
    idx_dict : dict
        Index conversion dictionary outputted by :code:`get_idx_dict`
    seq_dict : dict
        Maps the residue index to three letter amino acid
    sel_chain_id : str
        The chain to be designed. All residues outside of this chain will be treated as fixed to the native sequence.

    Returns
    =======
    subset_etab_matrix : np.ndarray
        The new etab
    subset_E_idx : np.ndarray
        The new indexing matrix
    """
    # etab will have dimensions of l x k x 20 x 20
    # l is the number of residues in selChainID
    # k is the number of nearest neighbors (30 by default)

    # create a new etab_matrix with proper dimensions
    filtered_idx_dict = {}
    renumbered_idx_dict = {}
    new_idx = {} #original idx to idx after selecting the chain of interest
    for idx,(chain, resid) in idx_dict.items():
        if (chain == sel_chain_id) or (not sel_chain_id):
            new_idx[idx] = len(filtered_idx_dict)
            filtered_idx_dict[idx] = (chain,resid)
            renumbered_idx_dict[new_idx[idx]] = (chain,resid)
    l = len(filtered_idx_dict)
    k = etab_matrix.shape[1] #normally 30
    subset_etab = np.zeros((l,k,20,20),dtype=float)
    # print('idx_dict',idx_dict)
    # print('filtered_idx_dict',filtered_idx_dict)
    # print('new_idx',new_idx)

    # create a new E_idx with the new residue indices
    subset_E_idx = np.zeros((l,k),dtype=float)
    for res_i_idx,neighbors in enumerate(E_idx):
        if (res_i_idx in filtered_idx_dict) or (not sel_chain_id):
            for nth_neighbor,res_j_idx in enumerate(neighbors):
                idx = new_idx[res_j_idx] if res_j_idx in new_idx else -1
                subset_E_idx[new_idx[res_i_idx],nth_neighbor] = idx

    # print('E_idx',E_idx)
    # print('subset_E_idx',subset_E_idx)

    # consider the following scenarios for a pair of residues (i,j)
    # 1) both are selected: copy all energies involving this pair of residues to the new etab
    # 2) i is selected and j is not:
    #   A) copy the self energies of i, discard the self energies of j. 
    #   B) take the pair energy of E(i,a,j,b), where a is one of 20 amino acids at position i and b is the native amino acid from the PDB file at position j, and add to the self energy E(i,a)
    # 3) both are not selected: omit all energies corresponding to these positions (as k is fixed to 30, this will mean leaving energies set to 0)
    # final note: since there are two pair energies for each E(i,a,j,b), then we must average them before adding to the self energies

    variable_fixed_pair = {} # dict[(var_res_idx,fixed_res_idx)] = energies : np.ndarray 20 x 20
    for i_idx,pos_i_slice in enumerate(etab_matrix):
        i_in_subset = (i_idx in filtered_idx_dict)
        if i_in_subset:
            # scenario 1 or 2A
            subset_etab[new_idx[i_idx]][0] = pos_i_slice[0]

        for nth_neighbor,ij_energies in enumerate(pos_i_slice):
            j_in_subset = (E_idx[i_idx][nth_neighbor] in filtered_idx_dict)
            if i_in_subset and j_in_subset:
                # scenario 1
                subset_etab[new_idx[i_idx]][nth_neighbor] = ij_energies

            elif i_in_subset or j_in_subset:
                # scenario 2B
                # the original residue idxs are unambiguously ordered based on whether they are selected: (variable,fixed)
                if i_in_subset: #i is variable
                    idxs = (i_idx,E_idx[i_idx][nth_neighbor])
                    if idxs in variable_fixed_pair:
                        variable_fixed_pair[idxs] = np.mean(np.array([variable_fixed_pair[idxs],ij_energies]),axis=0)
                    else:
                        variable_fixed_pair[idxs] = ij_energies
                else: #j is variable
                    idxs = (E_idx[i_idx][nth_neighbor],i_idx)
                    # before averaging, the energies are transposed so that they both have the "perspective" of the variable residue (variable residue should be i)
                    if idxs in variable_fixed_pair:
                        variable_fixed_pair[idxs] = np.mean(np.array([variable_fixed_pair[idxs],np.transpose(ij_energies)]),axis=0)
                    else:
                        variable_fixed_pair[idxs] = np.transpose(ij_energies)

            else:
                continue

    # scenario 2B: add pair to self
    for (var_idx,fix_idx),energies in variable_fixed_pair.items():
        # get the amino acid idx at the fixed position
        if seq_dict[fix_idx] in AA_to_int:
            fixed_aa_idx = AA_to_int[seq_dict[fix_idx]]
        else:
            residue = aa_to_similar(seq_dict[fix_idx])
            if residue not in AA_to_int:
                raise ValueError(residue+' not recognized')
            fixed_aa_idx = AA_to_int[residue]

        # extract the 20 pair energies corresponding to fixed_aa_idx and diagonalize
        self_e_corr = np.diagflat(energies[:,fixed_aa_idx])
        
        # if (var_idx == 35):
        #     print(var_idx,fix_idx)
        #     print(seq_dict[fix_idx])
        #     print(AA_to_int[seq_dict[fix_idx]])
        #     print(energies[:,fixed_aa_idx])
        
        # add energies to self at variable position
        subset_etab[new_idx[var_idx]][0] = subset_etab[new_idx[var_idx]][0] + self_e_corr

    return (subset_etab,subset_E_idx,renumbered_idx_dict)

# should work for multi-chain proteins now
def to_etab_file(etab_matrix, E_idx, idx_dict, out_path):
    """Write an :code:`.etab` file based on the fed in matrix and other indexing factors.

    Args
    ====
    etab_matrix : np.ndarray
        Etab outputted by TERMinator
    E_idx : np.ndarray
        Indexing matrix associated with :code:`etab_matrix`
    idx_dict : dict
        Index conversion dictionary outputted by :code:`get_idx_dict`
    out_path : str
        Path to write the etab to

    Returns
    =======
    bool
        Whether or not the parsing occured without errors
    out_path : str
        The output path fed in
    """
    out_file = open(out_path, 'w')

    # etab matrix: l x k x 20 x 20
    self_etab = etab_matrix[:, 0]
    pair_etab = etab_matrix[:, 1:]
    E_idx = E_idx[:, 1:]

    # l x 20
    self_nrgs = np.diagonal(self_etab, offset=0, axis1=-2, axis2=-1)
    for aa_idx, aa_nrgs in enumerate(self_nrgs):
        # pylint: disable=broad-except
        try:
            chain, resid = idx_dict[aa_idx]
        except Exception:
            eprint("num residues: ", len(self_nrgs))
            eprint(out_path)
            eprint(idx_dict)
            traceback.print_exc()
            return False, out_path
        for aa_int_id, nrg in enumerate(aa_nrgs):
            aa_3lt_id = int_to_3lt_AA[aa_int_id]
            out_file.write('{},{} {} {}\n'.format(chain, resid, aa_3lt_id, nrg))

    pair_nrgs = {}

    # l x k-1 x 20 x 20
    for i_idx, nrg_slice in enumerate(pair_etab):
        for k, k_slice in enumerate(nrg_slice):
            j_idx = E_idx[i_idx][k]
            if j_idx == -1:
                continue
            chain_i, i_resid = idx_dict[i_idx]
            chain_j, j_resid = idx_dict[j_idx]

            for i, i_slice in enumerate(k_slice):
                i_3lt_id = int_to_3lt_AA[i]
                for j, nrg in enumerate(i_slice):
                    j_3lt_id = int_to_3lt_AA[j]

                    # every etab has two entries i, j and j, i
                    # average these nrgs
                    key = [(chain_i, i_resid, i_3lt_id), (chain_j, j_resid, j_3lt_id)]
                    key.sort(key=lambda x: x[1])
                    key = tuple(key)
                    if key not in pair_nrgs.keys():
                        pair_nrgs[key] = nrg
                    else:
                        current_nrg = pair_nrgs[key]
                        pair_nrgs[key] = (current_nrg + nrg) / 2

    for key, nrg in sorted(pair_nrgs.items(), key=lambda pair: pair[0][0][1]):
        chain_i, i_resid, i_3lt_id = key[0]
        chain_j, j_resid, j_3lt_id = key[1]
        out_file.write('{},{} {},{} {} {} {}\n'.format(chain_i, i_resid, chain_j, j_resid, i_3lt_id, j_3lt_id, nrg))

    out_file.close()
    return True, out_path


def get_idx_dict(pdb, chain_filter=None):
    """From a :code:`.red.pdb` file, generate a dictionary mapping indices used within TERMinator
    to indices used by the :code:`.red.pdb` file.

    Args
    ====
    pdb : str
        path to :code:`.red.pdb` file
    chain_filter : list of str or None
        only parse chains from :code:`chain_filter`. If :code:`None`, parse
        all chains

    Returns
    =======
    idx_dict : dict
        Dictionary mapping indices used within TERMinator
        to indices used by the :code:`.red.pdb` file.
    seq_dict : dict
        Maps residue index to amino acid index
    """
    idx_dict = {}
    seq_dict = {}
    with open(pdb, 'r') as fp:
        current_idx = 0
        for line in fp:
            data = line.strip()
            if data == 'TER' or data == 'END':
                continue
            try:
                chain = data[21]
                # residx = int(data[22:26].strip())
                # icode = data[26]
                # if icode != ' ':
                #     residx = str(residx) + icode
                residx = data[22:27].strip()  # rip i didn't know about icodes
                resname = data[17:20]

            except Exception as e:
                print(data)
                raise e

            if chain_filter:
                if chain not in chain_filter:
                    continue

            if (chain, residx) not in idx_dict.values():
                idx_dict[current_idx] = (chain, residx)
                seq_dict[current_idx] = resname
                current_idx += 1

    return (idx_dict,seq_dict)


def _to_numpy_file_wrapper(filepath, out_folder):
    """Wrapper for :code:`parseEtab` for path manipuation and error catching.
    Args
    ----
    in_path : str
        input .etab to :code:`parseEtab`
    out_folder : str
        output directory to dump .etab.npy files into
    """
    name = re.split('/|\\\\', filepath)[-1]
    out_file = os.path.join(out_folder, name)
    print('out file', out_file)
    try:
        _, _, etab = parseEtab(filepath, save=False)
        np.save(out_file, etab)
    except Exception as e:
        print(out_file, file=sys.stderr)
        print(e)
        raise e
    return True, out_file


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Generate etabs')
    parser.add_argument('--output_dir', help='output directory', required=True)
    parser.add_argument("--dtermen_data", help="Root directory for all dTERMen runs", required=False)
    parser.add_argument('--num_cores', help='number of processes for parallelization', default=1)
    parser.add_argument('--chain_id', help='if provided, will selected a subset of the etab positions with this chain ID as variable and set the rest to fixed',default='')
    parser.add_argument('--chain_id_from_name', help='if provided, will extract the chain ID from the name of the structure (e.g. if name is PDBIB_A_B_1_2, chain ID will be B)',action='store_true')
    parser.add_argument('--save_type', help='file format to save in', default='etab')
    parser.add_argument('-u', dest='update', help='flag for force updating etabs', default=False, action='store_true')
    args = parser.parse_args()
    args.update = True # JFM
    if not os.path.isdir(os.path.join(args.output_dir, 'etabs')):
        os.mkdir(os.path.join(args.output_dir, 'etabs'))

    print(f"results dir: {args.output_dir}")

    if args.save_type == "etab":
        with open(os.path.join(args.output_dir, 'net.out'), 'rb') as fp:
            dump = pickle.load(fp)
        print("loaded dump")
        pbar = tqdm(total=len(dump))
    else:
        pbar = tqdm(total=47) ## JFM

    pool = mp.Pool(int(args.num_cores))
    start = time.time()
    not_worked = []

    def check_worked(res):
        """Update progress bar per iteration"""
        worked, out_path = res
        pbar.update()
        if not worked:
            not_worked.append(out_path)

    def raise_error(error):
        """Propogate error upwards"""
        raise error
    if args.save_type == "etab":
        print("starting etab dump")
        for data in dump:
            pdb = data['ids'][0]
            E_idx = data['idx'][0].copy()
            etab = data['out'][0].copy()

            chain_id = args.chain_id
            if args.chain_id_from_name:
                chain_id = pdb.split('_')[2]

            print(pdb,chain_id)
            pdb_path = find_dtermen_folder(pdb, args.dtermen_data)
            (idx_dict,seq_dict) = get_idx_dict(os.path.join(pdb_path, f'{pdb}.red.pdb'))
            (subset_etab,subset_E_idx,renumbered_idx_dict) = select_variable_positions(etab,E_idx,idx_dict,seq_dict,chain_id)

            out_path = os.path.join(args.output_dir, 'etabs/' + pdb + '.etab')

            if os.path.exists(out_path) and not args.update:
                print(f"{pdb} already exists, skipping")
                pbar.update()
                continue


            res = pool.apply_async(_to_etab_file_wrapper,
                                args=(subset_etab, subset_E_idx, renumbered_idx_dict, out_path),
                                callback=check_worked,
                                error_callback=raise_error)
        else:
            args.output_dir = os.path.join(args.output_dir, "etabs")
            in_files = glob.glob(args.output_dir+"/*.etab")
            in_list = [os.path.abspath(i) for i in in_files]
            for filepath in in_list:
                res = pool.apply_async(_to_numpy_file_wrapper,
                                    args=(filepath, args.output_dir),
                                    callback=check_worked,
                                    error_callback=raise_error)

    pool.close()
    pool.join()
    pbar.close()
    print(f"errors in {not_worked}")
    for path in not_worked:
        os.remove(path)
    end = time.time()
    print(f"done, took {end - start} seconds")
