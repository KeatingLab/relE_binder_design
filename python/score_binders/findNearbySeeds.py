import argparse
import json
from multiprocessing.sharedctypes import Value
import os,io,time
import pickle
import numpy as np
import sys

import torch
import torch.nn as nn
from torch.utils.data import DataLoader

from scripts.data.preprocessing.cleanStructs import extractBackbone
from terminator.data.data import BinderScoringIterableDataset,ComplexScoringDataset
# from terminator.models.TERMinator import TERMinator
# from terminator.utils.model.loop_utils import run_epoch,_to_dev
# from terminator.utils.model.loss_fn import construct_loss_fn
# from terminator.utils.common import int_to_3lt_AA,AA_to_int,int_to_AA
# from terminator.utils.model.default_hparams import DEFAULT_MODEL_HPARAMS, DEFAULT_TRAIN_HPARAMS

# from python.score_binders.score_utils import *
from python.utils.utils import *

def load_hparams(model_path, default_hparams, output_name):
    print("loading params")
    # load hparams
    hparams_path = os.path.join(model_path, output_name)
    hparams = json.load(open(hparams_path, 'r'))
    for hparam in default_hparams:
        if hparam not in hparams:
            # print(f"{hparam} = {default_hparams[hparam]}")
            hparams[hparam] = default_hparams[hparam]
    return hparams

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Score peptide-binder backbones using TERMinator energies')
    parser.add_argument('--binder_dataset', help='Multi-entry PDB file containing a target structure and various binder structures')
    # parser.add_argument('--binder_list', help='A file where each line the name of a binder in the binder_dataset file that should be scored', default='')
    parser.add_argument('--pdb',help='pdb file describing the structure of the target protein', default='')
    parser.add_argument('--sel_chain_ids', help='the chain IDs to select from the structure, ex: "--sel_chain_ids A,B"', required=True, type=str)

    args = parser.parse_args()
    print(sys.argv)

    # Initialize the special dataloader for binders
    dataset = BinderScoringIterableDataset(args.binder_dataset,args.pdb,1000,None,mode="binder_only",skip_package=True)

    sel_chain_ids = args.sel_chain_ids.split(',')
    print(f"selected chain IDs {sel_chain_ids}")

    findSeedsCloseToSelection(args.pdb,sel_chain_ids,dataset,distance_cutoff=3.5)

    print('Done!')
