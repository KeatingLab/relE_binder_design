"""Run NRGTEN to calculate dynamic signatures for all proteins in input list.

Usage:
    .. code-block::

        python generateNrgtenData.py \\
            --pdb_root \\
            --pdb_list \\
            --out_dir \\
            --nrgten_type \\
            --num_processes \\
            --success_file \\

"""

import argparse
import os
import multiprocessing as mp
from nrgten.encom import ENCoM
import matplotlib.pyplot as plt
from tqdm import tqdm
import pickle
import numpy as np
from parseCoords import parseCoords

def nrgten_calc(pdb, nrgten_type, out_dir, update):
    """Calculates nrgten features for given pdb file

    Args
    ----
    pdb : str
        ID of current pdb file to load
    nrgten_type : str
        Type of nrgten calculation to perform
    out_dir : str
        Path to directory containing results
    update : bool
        Whether to update existing files

    Returns
    -------
    nrgten_success : bool
        Indicates whether nrgten calculation was successful
    pdb : str
        ID of current pdb file to load
    """
    pdb = pdb.strip()
    pdb_file = os.path.join(pdb, pdb + ".red.pdb")
    if not os.path.exists(pdb_file):
        print(f"No file {pdb_file} exists.")
        return False, None
    try:
        print(pdb_file)
        model = None
        if nrgten_type.find("signatures") > -1:
            nrgten_path = os.path.join(out_dir, pdb[:-1] + "_dynamic_signatures.txt")
            if os.path.exists(nrgten_path) and not update:
                if nrgten_type.find("ensembles") == -1:
                    return True, pdb
            else:
                model = ENCoM(pdb_file)
                dynamic_sigs = model.compute_bfactors()
                with open(nrgten_path, 'w') as of:
                    for sig in dynamic_sigs:
                        of.write(str(sig) + "\n")
                of.close()
        if nrgten_type.find("ensembles") > -1:
            nrgten_path = os.path.join(out_dir, pdb + "_conformational_ensembles.pdb")
            nrgten_features_path = os.path.join(out_dir, pdb + "_conformational_ensembles.features")
            if os.path.exists(nrgten_path) and os.path.exists(nrgten_features_path) and not update:
                return True, pdb
            if model is None:
                model = ENCoM(pdb_file)
            if not os.path.exists(nrgten_path) or update:
                model.build_conf_ensemble([7, 8], nrgten_path, step=2, max_displacement=2, )
            coords, _ = parseCoords(nrgten_path, ensemble=True, save=False, verbose=False)
            coords_tensor = None
            if len(coords) == 1:
                chain = next(iter(coords.keys()))
                coords_tensor = coords[chain]
            else:
                chains = sorted(coords.keys())
                coords_tensor = np.vstack([coords[c] for c in chains])
            output = {
                'coords': coords_tensor,
            }
            with open(nrgten_features_path, 'wb') as fp:
                pickle.dump(output, fp)
            print("made it to end!")

        return True, pdb
    except Exception as e:
        print(f"Error {e} for pdb {pdb}.")
        return False, pdb


if __name__ == '__main__':
    # TODO: create argument for parser to take in folder w proper dir structure for pdbs
    print("starting main")
    parser = argparse.ArgumentParser('Run NRGTEN to calculate dynamic signatures.')
    parser.add_argument("--pdb_root", help="The root for all raw data PDB files", required=True)
    parser.add_argument("--pdb_list", help="The list of pdb files to use", required=True)
    parser.add_argument('--out_dir', help='Output directory', required=True)
    parser.add_argument('--nrgten_type', help="Type of nrgten calculation to perform", required=True)
    parser.add_argument('--num_processes', help="Number of processes for parallelism", required=True)
    parser.add_argument('--success_file', help="File to save successful pdbs", required=True)
    parser.add_argument('--update',help='Whether to update existing files', default=False)
    args = parser.parse_args()
    print(args.pdb_root)

    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir)
    os.chdir(args.pdb_root)

    pdbs = []
    with open(args.pdb_list, 'r') as f:
        pdbs = f.readlines()
    print(f"starting pool for {len(pdbs)} pdbs.")
    nrgten_successes = []
    with mp.Pool(int(args.num_processes)) as pool:
        progress = tqdm(total=len(pdbs))

        def update_progress(res):
            del res
            progress.update(1)

        res_list = [
            pool.apply_async(nrgten_calc, (id, args.nrgten_type, args.out_dir, args.update),
                                callback=update_progress) for id in pdbs
        ]
        pool.close()
        pool.join()
        progress.close()

        for res in res_list:
            data = res.get()
            if data is not None and data[0]:
                nrgten_successes.append(data[1])
            elif data and not data[0]:
                print(data[1])


    outfile = os.path.join(args.out_dir, args.success_file)
    with open(outfile, 'w') as of:
        for pdb in nrgten_successes:
            of.write(pdb)
    of.close()
