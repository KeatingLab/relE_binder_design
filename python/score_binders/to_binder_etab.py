import argparse
import json
import os,io,time
import pickle
import numpy as np

import torch
import torch.nn as nn
from torch.utils.data import DataLoader

from scripts.data.preprocessing.cleanStructs import extractBackbone
from terminator.data.data import BinderScoringIterableDataset,ComplexScoringDataset
from terminator.models.TERMinator import TERMinator
from terminator.utils.model.loop_utils import run_epoch,_to_dev
from terminator.utils.model.loss_fn import construct_loss_fn
from terminator.utils.common import int_to_3lt_AA,AA_to_int,int_to_AA

def reduce_etab_to_variable_positions(etab, E_idx, seq, peptide_res):
    peptide_range = (peptide_res.min(),peptide_res.max()) # p = peptide residue length

    # get the subset of the energy table that is centered around peptide residues 
    variable_etab = etab[peptide_range[0]:peptide_range[1]] # p x k
    variable_E_idx = E_idx[peptide_range[0]:peptide_range[1]] # p x k

    # for every interaction between a peptide residue and a protein residue, get the energy given the protein residue native aa identity and add to the self energy
    variable_etab_selfE = torch.diagonal(variable_etab[:,0].squeeze(1),offset=0,dim1=-2,dim2=-1) # p x 20
    
    variable_etab_pairE = variable_etab[:,1:] # p x 29 x 20 x 20
    variable_E_idx_pair = E_idx[:,1:] # p X 29

    # mask the pairs according to whether the interacting residue is protein or peptide
    mask_unfixed = np.isin(variable_E_idx_pair, peptide_res) # p x 29

    variable_etab_pairE_fixed = variable_E_idx_pair * ~mask_unfixed
    variable_etab_pairE_unfixed = variable_E_idx_pair * mask_unfixed

    variable_E_idx_pair_expanded = np.expand_dims(variable_E_idx_pair,(2,3)). # p x 29 x 20 x 20
    variable_etab_pairE_fixed_u = np.take(variable_etab_pairE_fixed)


    return variable_etab, variable_E_idx

def write_etab_to_file():
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Load network output from scoreBinders.py and convert to dTERMen format energy table')
    parser.add_argument('--complex_netout_list', help='A file where each line is the path to a network output pickle file')
    args = parser.parse_args()

    complex_netout_paths = []
    with open(args.complex_netout_list,'r') as file:
        complex_netout_paths = [line.rstrip() for line in file]
    
    for path in complex_netout_paths:
        with open(path,'rb') as file:
            complex_netout = pickle.load(file)

        



    print('Done!')