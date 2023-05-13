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

## for loading structures from multi-entry PDB file and threading a new sequence on
def getPDBsFromMultiEntryFile(pdbNames,multientryPDBpath,outDir):
    pdbNames_set = set(pdbNames)
    pdbNames_found = set()
    print(f"Searching for {len(pdbNames_set)} structures in multi-entry PDB file")
    # open multi-entry pdb file
    with open(multientryPDBpath,'r') as file:
        # get lines corresponding to selected PDBs
        filenames = []
        for line in file:
            if line[:6] == 'HEADER':
                name = line[10:].rstrip()
                if name not in pdbNames_set:
                    continue
                if name in pdbNames_found:
                    print(f"Warning: {name} was observed twice, there are duplicate names in the PDB file")
                pdbNames_found.add(name)
                pdb_lines = []
                pdb_line = file.readline()
                while pdb_line[:3] != 'END':
                    if pdb_line[:3] != 'TER':
                        pdb_lines.append(pdb_line)
                    pdb_line = file.readline()
                with open(os.path.join(outDir,name+'.pdb'),'w') as pdb_file:
                    for pdb_line in pdb_lines:
                        pdb_file.write(pdb_line)
                filenames.append(name)

    print(f"Found {len(pdbNames_found)} structures in the file")
    # return the list of paths
    return filenames
    

def addPDBsToSilentFile(filtered_df,multientryPDBpath,silent_binary_name,
                        targetPDBchains,batch_size,
                        binder_chain_id='0',name_col='name',seq_col='',
                        silentfrompdbs_path='/home/gridsan/sswanson/local_code_mirror/silent_tools/silentfrompdbs'):
    '''
    Load PDBs from multi-entry file, replace info, and add to silent binary file

    NOTE: function is optimized to limit the number of PDB files in a given directory 
    '''
    warnings.filterwarnings('ignore')
    io = PDBIO()
    p = PDBParser(PERMISSIVE=1)

    pdb_names_batch = filtered_df[name_col].unique()
    # extract the PDBs of interest from the multi-entry file
    originalPDBDir = silent_binary_name + '_temp/extractedPDBs'
    shutil.rmtree(originalPDBDir,ignore_errors=True)
    pathlib.Path(originalPDBDir).mkdir(parents=True, exist_ok=True)
    filenames = getPDBsFromMultiEntryFile(pdb_names_batch,multientryPDBpath,originalPDBDir)

    # load with Bio.PDB parser
    modifiedPDBDir = silent_binary_name + '_temp/modifiedPDBs'
    shutil.rmtree(modifiedPDBDir,ignore_errors=True)
    pathlib.Path(modifiedPDBDir).mkdir(parents=True, exist_ok=True)
    warnings.filterwarnings('ignore')
    for i,row in filtered_df[filtered_df[name_col].isin(pdb_names_batch)].iterrows():
        name = row[name_col]
        new_name = row[name_col] + '_' + row['unique_id']
        binder_structure = p.get_structure(name,os.path.join(originalPDBDir,name+'.pdb'))
        
        # set sequence/bfactor
        if seq_col != '':
            for i,residue in enumerate(binder_structure[0][binder_chain_id]):
                residue.resname = one_to_three(row[seq_col][i])
                
        # combine with target chain 
        if targetPDBchains:
            for i,chain in enumerate(targetPDBchains):
                binder_structure[0].add(chain)
        
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
    parser.add_argument('--multientry_pdb', type=str, help='')
    parser.add_argument('--silent_name', type=str, help='')
    parser.add_argument('--target_pdb', type=str, help='')
    parser.add_argument('--n_batches', type=int, help='')
    parser.add_argument('--batch_id', type=int, help='')
    parser.add_argument('--binder_chain_id', type=str, help='')
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

    # load target pdb
    p = PDBParser(PERMISSIVE=1)
    target_structure = p.get_structure('',args.target_pdb)
    targetPDBchains = [x for x in target_structure[0].get_chains()]

    silent_name = args.silent_name + "_" + str(args.batch_id)

    addPDBsToSilentFile(batch_df,args.multientry_pdb,silent_name,
                        targetPDBchains,batch_size,
                        binder_chain_id=args.binder_chain_id,name_col=args.name_col,seq_col=args.seq_col,
                        silentfrompdbs_path=args.silentfrompdbs)

    print("Done!")