"""Functions to parse :code:`.red.pdb` files"""
import pickle
from pyexpat import model

import numpy as np
import copy

from terminator.utils.common import aa_three_to_one,aa_to_similar
    
def parseCoords(filename, ensemble=False, save=True, verbose=True, record_b_factor=False, valid_entry_lines=None):
    """ Parse coordinates from :code:`.red.pdb` files  or a list of lines from extractBackbone, and dump in
    files if specified.

    Args
    ====
    filename : str
        path to :code:`.red.pdb` file
    ensemble : boot, default=False
        whether or not pdb file contains multiple conformations
    save : bool, default=True
        whether or not to dump the results
    verbose : bool, default=True
        whether or not to print intermediate results
    record_b_factor : bool, default=False
        whether or not to return B-factors

    valid_entry_lines : list, default=None
        PDB file lines from extractBackbone()

    Returns
    =======
    chain_tensors : dict
        Dictionary mapping chain IDs to arrays of atomic coordinates.

    seq : str
        Sequence of all chains concatenated.
    """
    chain_dict = {}
    chain_dict_save = {}
    if record_b_factor:
        b_factor_dict = {}
        b_factor_dict_save = {}
    if ensemble:
        model_dict = {}
    if filename is not None:
        with open(filename, 'r') as fp:
            lines = fp.readlines()
    elif valid_entry_lines is not None:
        lines = valid_entry_lines
    else:
        raise ValueError('Must provide either filename or valid entry lines')
    
    for line in lines:
        data = line.strip()
        if data[:3] == 'TER' or (data[:3] == 'END' and data[3:6] != 'MDL') or data[:6] == 'REMARK':
            continue
        elif data[:5] == 'MODEL':
            cur_model = data[11:14].strip()
            continue
        elif data[:6] == 'ENDMDL':
            if ensemble:
                model_dict[cur_model] = chain_dict
            if not chain_dict_save:
                chain_dict_save = copy.deepcopy(chain_dict)
            chain_dict = {}
            if record_b_factor:
                if not b_factor_dict_save:
                    b_factor_dict_save = copy.deepcopy(b_factor_dict)
                b_factor_dict = {}
            continue
        try:
            element = data[13:16].strip()
            residue = data[17:20].strip()
            residx = data[22:27].strip()
            chain = data[21]
            x = data[30:38].strip()
            y = data[38:46].strip()
            z = data[46:54].strip()
            coords = [float(coord) for coord in [x, y, z]]
            if record_b_factor:
                b_factor = float(data[61:67])
        except Exception as e:
            print(data)
            raise e

        if chain not in chain_dict.keys():
            chain_dict[chain] = {element: [] for element in ['N', 'CA', 'C', 'O']}
            chain_dict[chain]["seq_dict"] = {}
            chain_dict[chain]["res_list"] = []
            if record_b_factor:
                b_factor_dict[chain] = {element: [] for element in ['N', 'CA', 'C', 'O']}


        # naively model terminal carboxylate as a single O atom
        # (i cant find the two oxygens so im just gonna use OXT)
        if element == 'OXT':
            element = 'O'
        chain_dict[chain][element].append(coords)
        if record_b_factor:
            b_factor_dict[chain][element].append(b_factor)

        seq_dict = chain_dict[chain]["seq_dict"]
        res_list = chain_dict[chain]["res_list"]
        if residx not in seq_dict.keys():
            residue = aa_to_similar(residue)
            if verbose and (residue == '0Q4') or (residue == 'FOL'):
                print(filename, residue)
                
            seq_dict[residx] = aa_three_to_one(residue)
            res_list+=[(chain,residx)]
    chain_tensors = {}
    if record_b_factor:
        b_factor_tensors = {}
    seq = ""
    res_info = list()
    chain_ids = sorted(chain_dict.keys())
    if ensemble:
        chain_dict = chain_dict_save
    for chain in chain_ids:
        if ensemble:
            coords = [ [model_dict[model][chain][element] for model in model_dict.keys()] for element in ['N', 'CA', 'C', 'O'] ]
        else:
            coords = [chain_dict[chain][element] for element in ['N', 'CA', 'C', 'O']]
            if record_b_factor:
                b_factors = [b_factor_dict[chain][element] for element in ['N', 'CA', 'C', 'O']]
        chain_tensors[chain] = np.stack(coords, 1) if not ensemble else np.moveaxis(np.swapaxes(np.array(coords), 0, 2), 1, 3)
        if record_b_factor:
            b_factor_tensors[chain] = np.stack(b_factors, 1)
        seq_dict = chain_dict[chain]["seq_dict"]        
        chain_seq = "".join([seq_dict[i] for i in seq_dict.keys()])
        if verbose:
            print('chain seq: ',chain_seq)
        res_info += chain_dict[chain]["res_list"]
        assert len(chain_seq) == chain_tensors[chain].shape[0], (chain_seq, chain_tensors[chain].shape, filename)
        seq += "".join([seq_dict[i] for i in seq_dict.keys()])

    if save:
        with open(filename[:-8] + '.coords', 'wb') as fp:
            pickle.dump(chain_tensors, fp)
        with open(filename[:-8] + '.seq', 'w') as fp:
            fp.write(seq)
        if record_b_factor:
            with open(filename[:-8] + '.bfactors', 'wb') as fp:
                pickle.dump(b_factor_tensors, fp)

    if not record_b_factor:
        return chain_tensors, seq, res_info
    return b_factor_tensors
