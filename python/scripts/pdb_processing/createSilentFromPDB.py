import argparse
import os
import glob
import json
import pickle
import pathlib
import subprocess
import shutil
import warnings
import sys

import pandas as pd
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import one_to_three 
from Bio.PDB import PDBIO

import numpy as np
import pandas as pd

def PDBtoSilentFile(filtered_df,
                    PDB_dir,
                    silent_binary_name,
                    name_col='name',
                    seq_col='mcmc_seq_peptide',
                    binder_chain_id='0',
                    silentfrompdbs_path='/home/gridsan/sswanson/local_code_mirror/silent_tools/silentfrompdbs'):
    
    warnings.filterwarnings('ignore')
    io = PDBIO()
    p = PDBParser(PERMISSIVE=1)

    modifiedPDBDir = silent_binary_name + '_temp/modifiedPDBs'
    shutil.rmtree(modifiedPDBDir,ignore_errors=True)
    pathlib.Path(modifiedPDBDir).mkdir(parents=True, exist_ok=True)

    pdb_names_batch = filtered_df[name_col].unique()
    for i,row in filtered_df[filtered_df[name_col].isin(pdb_names_batch)].iterrows():
        name = row[name_col]
        new_name = row[name_col] + '_' + row['unique_id']
        binder_structure = p.get_structure(name,os.path.join(PDB_dir,name+'.pdb'))
                
        # # combine with target chain 
        # if targetPDBchains:
        #     for i,chain in enumerate(targetPDBchains):
        #         binder_structure[0].add(chain)

        for i,residue in enumerate(binder_structure[0][binder_chain_id]):
            residue.resname = one_to_three(row[seq_col][i])
        
        # write modified file
        io.set_structure(binder_structure)
        io.save(os.path.join(modifiedPDBDir,new_name+'.pdb'))

    # external call to load modified PDBs to silent file
    command = f"{silentfrompdbs_path} {modifiedPDBDir}/*.pdb > {silent_binary_name}_binder_designs.silent"
    CompletedProcess = subprocess.run(command,shell=True)
    if CompletedProcess.returncode != 0:
        raise ValueError('Subprocess returned an error when trying to create a silent file')

    # delete PDBs and work on next batch
    shutil.rmtree(silent_binary_name + '_temp')

    warnings.filterwarnings('default')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Load sequence and peptide structures and create rosetta silent files')
    parser.add_argument('--selected_designs', type=str, help='')
    parser.add_argument('--pdb_dir', type=str, help='')
    parser.add_argument('--silent_name', type=str, help='')
    parser.add_argument('--n_batches', type=int, help='')
    parser.add_argument('--batch_id', type=int, help='')
    parser.add_argument('--name_col', type=str, help='')
    parser.add_argument('--seq_col', type=str, help='')
    parser.add_argument('--silentfrompdbs', type=str, help='', default='/home/gridsan/sswanson/local_code_mirror/silent_tools/silentfrompdbs')
    args = parser.parse_args()
    print(sys.argv)

    # get section of df
    sel_df = pd.read_csv(args.selected_designs)
    sel_df['unique_id'] = sel_df.index
    sel_df['unique_id'] = sel_df['unique_id'].astype(str)
    batch_size = len(sel_df) // (args.n_batches-1)
    batch_df = sel_df.iloc[args.batch_id*batch_size:args.batch_id*batch_size+batch_size]

    # # load target pdb
    # p = PDBParser(PERMISSIVE=1)
    # target_structure = p.get_structure('',args.target_pdb)
    # targetPDBchains = [x for x in target_structure[0].get_chains()]

    silent_name = args.silent_name + "_" + str(args.batch_id)

    PDBtoSilentFile(batch_df,
                    args.pdb_dir,
                    silent_name,
                    name_col=args.name_col,
                    seq_col=args.seq_col,
                    binder_chain_id='0',
                    silentfrompdbs_path=args.silentfrompdbs)

    print("Done!")