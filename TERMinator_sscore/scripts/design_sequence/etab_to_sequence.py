# load TERMinator dependencies

# load PuLP/CPLEX dependencies

"""Use ILP to optimize an energy table and generate a sequence for a given backbone structure

Usage:
    .. code-block::

        python etab_to_sequence.py \\
            # --netout <dataset_dir> \\
            # [--charge_constraint <trained_model_dir>] \\

See :code:`python etab_to_sequence.py --help` for more info.
"""

import argparse
import json
from multiprocessing.sharedctypes import Value
import os,io,time
import pickle
import numpy as np
import sys
from tqdm import tqdm
import multiprocessing as mp
import traceback

import torch

from terminator.utils.common import int_to_3lt_AA,AA_to_int,int_to_AA,aa_three_to_one
from scripts.data.postprocessing.search_utils import find_dtermen_folder
from scripts.data.postprocessing.to_etab import get_idx_dict,eprint
from scripts.design_sequence.ilp_utilities import optimizeEnergyTableViaLP

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Use ILP to optimize an energy table and generate a sequence design')
    parser.add_argument('--output_dir', help='The directory with netout from TERMinator', required=False, default="")
    parser.add_argument("--dtermen_data", help="Root directory for all dTERMen runs", required=False, default="")
    parser.add_argument("--dtermen_etabs", help="A file with paths to '*.etab' dTERMen-formatted energy tables", required=False, default="")
    parser.add_argument('--num_cores', help='number of processes for parallelization', default=1)
    # parser.add_argument('--chain_id', help='if provided, will selected a subset of the etab positions with this chain ID as variable and set the rest to fixed',default='')
    # parser.add_argument('--charge_constraint', help='The range of tolerable charges, given as [a,b] where a is most negative charge and b is the most positive charge tolerated, e.g., [-4,4]. Assumes physiological pH when determing charge state of each amino acid type.', action=argparse.store)
    args = parser.parse_args()
    args.update = True # JFM
    print(sys.argv)

    if (args.output_dir == "" and args.dtermen_etabs == ""):
        raise ValueError("Must provide '--output_dir' or '--dtermen_etabs'")

    # # parser = argparse.ArgumentParser('Generate etabs')
    # parser.add_argument('--output_dir', help='output directory', required=True)
    # parser.add_argument("--dtermen_data", help="Root directory for all dTERMen runs", required=False)
    # parser.add_argument('--num_cores', help='number of processes for parallelization', default=1)
    # parser.add_argument('--chain_id', help='if provided, will selected a subset of the etab positions with this chain ID as variable and set the rest to fixed',default='')
    # parser.add_argument('--chain_id_from_name', help='if provided, will extract the chain ID from the name of the structure (e.g. if name is PDBIB_A_B_1_2, chain ID will be B)',action='store_true')
    # parser.add_argument('--save_type', help='file format to save in', default='etab')
    # parser.add_argument('-u', dest='update', help='flag for force updating etabs', default=False, action='store_true')
    # args = parser.parse_args()
    # print(sys.argv)
    # args.update = True # JFM

    if args.output_dir != "":
        with open(os.path.join(args.output_dir, 'net.out'), 'rb') as fp:
            dump = pickle.load(fp)
        print("loaded netout")
    else:
        with open(args.dtermen_etabs,'r') as file:
            dump = [path.rstrip() for path in file]
        print("loaded dtermen etab paths")
    # pbar = tqdm(total=len(dump))

    # pool = mp.Pool(int(args.num_cores))
    start = time.time()
    # not_worked = []

    # def check_worked(res):
    #     """Update progress bar per iteration"""
    #     worked, out_path = res
    #     pbar.update()
    #     if not worked:
    #         not_worked.append(out_path)

    # def raise_error(error):
    #     """Propogate error upwards"""
    #     raise error

    opt = optimizeEnergyTableViaLP()

    print("Starting sequence design")
    for data in dump:
        if args.output_dir:
            pdb = data['ids'][0]
            E_idx = data['idx'][0]
            etab = data['out'][0]

            print(pdb)
            pdb_path = find_dtermen_folder(pdb, args.dtermen_data)
            (idx_dict,seq_dict) = get_idx_dict(os.path.join(pdb_path, f'{pdb}.red.pdb'))
            res_info = [(chain_id,resnum) for chain_id,resnum in idx_dict.values()]
            native_sequence = ''.join([aa_three_to_one(aa3) for aa3 in seq_dict.values()])

            opt.knn_etab_to_sequence(pdb, etab, E_idx, res_info, '', None, native_sequence)

            # res = pool.apply_async(_etab_to_sequence_wrapper,
            #                     args=(etab, E_idx, res_info, out_file, native_sequence),
            #                     callback=check_worked,
            #                     error_callback=raise_error)
        else:
            print("dTERMen etab: ",data)
            opt.dtermen_etab_to_sequence(data)

    # pool.close()
    # pool.join()
    # pbar.close()
    # print(f"errors in {not_worked}")
    # for path in not_worked:
    #     os.remove(path)
    end = time.time()
    print(f"done, took {end - start} seconds")