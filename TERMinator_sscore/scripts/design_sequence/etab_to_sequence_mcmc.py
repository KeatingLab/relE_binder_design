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
import pandas as pd
import sys
# from tqdm import tqdm
# import multiprocessing as mp
# import traceback

import torch

from terminator.utils.common import int_to_3lt_AA,AA_to_int,int_to_AA,aa_three_to_one
from scripts.data.postprocessing.search_utils import find_dtermen_folder
from scripts.data.postprocessing.to_etab import get_idx_dict,eprint
from scripts.design_sequence.mcmc_utilities import simulatedAnnealingMCMC, mcmcParams, residueGroup

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Use simulated annealing markov chain monte carlo to optimize an energy table and generate a sequence design')
    parser.add_argument('--output_dir', help='The directory with netout from TERMinator', required=False)
    parser.add_argument("--dtermen_data", help="Root directory for all dTERMen runs", required=False, default="")
    # parser.add_argument("--dtermen_etabs", help="A file with paths to '*.etab' dTERMen-formatted energy tables", required=False, default="")
    # parser.add_argument('--num_cores', help='number of processes for parallelization', default=1)
    parser.add_argument('--chain_id', help='if provided, will selected the etab positions with this chain ID as variable and set the rest to fixed',default='')
    parser.add_argument('--n_cyc', help='The number of independent mcmc cycles',default=100)
    parser.add_argument('--n_it', help='The number of iterations per mcmc cycle',default=1e6)
    parser.add_argument('--Ti',help='The initial temperature',default=1.0)
    parser.add_argument('--Tf',help='The final temperature',default=0.01)
    parser.add_argument('--seed',help='For reproducibility',default=42)
    parser.add_argument('--complexity_weight',help='If > 0, the complexity of the sequence will be calculated and this will be added to the energy to form a weighted sum',default=-1.0)
    parser.add_argument('--early_stopping',help='If > 0, will terminate early if the best sequence is repeatedly the same this many times',default=-1)
    args = parser.parse_args()
    args.update = True
    print(sys.argv)

    with open(os.path.join(args.output_dir, 'net.out'), 'rb') as fp:
        dump = pickle.load(fp)
    print("loaded netout")

    start = time.time()

    opt = simulatedAnnealingMCMC()
    params = mcmcParams(args.n_it,args.n_cyc,args.Ti,args.Tf,args.seed,args.complexity_weight,args.early_stopping)

    print("Starting sequence design")
    df_list = []
    for data in dump:
        # Load etab and residue info
        pdb = data['ids'][0]
        E_idx = data['idx'][0]
        etab = data['out'][0]
        print(pdb)
        pdb_path = find_dtermen_folder(pdb, args.dtermen_data)
        (idx_dict,seq_dict) = get_idx_dict(os.path.join(pdb_path, f'{pdb}.red.pdb'))
        res_info = [(chain_id,resnum) for chain_id,resnum in idx_dict.values()]
        native_sequence = ''.join([aa_three_to_one(aa3) for aa3 in seq_dict.values()])
        
        opt.loadEnergyTable(etab,E_idx,res_info)
        if args.chain_id:
            selResList = [residueGroup(args.chain_id)]
            opt.setConstraints(native_sequence,selResList)
        
        opt.sample(params)
        df = opt.getDataFrame()
        df['name'] = pdb
        df_list.append(df)
        
    df = pd.concat(df_list)
    df.to_csv('mcmc_optimization_results.csv')

    end = time.time()
    print(f"done, took {end - start} seconds")