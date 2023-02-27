"""Score protein backbones using TERMinator energy tables
"""

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
from terminator.data.data import SingleChainScoringDataset
from terminator.models.TERMinator import TERMinator
from terminator.utils.model.loop_utils import run_epoch,_to_dev
from terminator.utils.model.loss_fn import construct_loss_fn
from terminator.utils.common import int_to_3lt_AA,AA_to_int,int_to_AA
from terminator.utils.model.default_hparams import DEFAULT_MODEL_HPARAMS, DEFAULT_TRAIN_HPARAMS

from python.score_binders.score_utils_singlechain import *

# pylint: disable=unspecified-encoding

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
    parser = argparse.ArgumentParser('Score protein backbones using TERMinator energies and structure score')
    parser.add_argument('--structure_dataset', help='A file where each line is the path to a PDB file.')
    parser.add_argument('--model_dir', help='trained model folder', required=True)
    parser.add_argument('--seq_mode', help='controls which sequence is used when scoring the binder', default='consensus_aa')
    parser.add_argument('--dev', help='device to use', default='cuda:0')
    parser.add_argument('--store', help='if given, will write the output to a JSON file', action=argparse.BooleanOptionalAction)
    args = parser.parse_args()
    print(sys.argv)

    dev = args.dev
    if torch.cuda.device_count() == 0:
        dev = "cpu"

    # Initialize the special dataloader for binders
    dataset = SingleChainScoringDataset(args.structure_dataset)
    dataset_iter = iter(dataset)

    # Load the model configuration parameters
    model_hparams = load_hparams(args.model_dir, DEFAULT_MODEL_HPARAMS, 'model_hparams.json')
    run_hparams = load_hparams(args.model_dir, DEFAULT_TRAIN_HPARAMS, 'run_hparams.json')

    # backwards compatibility
    if "cov_features" not in model_hparams.keys():
        model_hparams["cov_features"] = False
    if "term_use_mpnn" not in model_hparams.keys():
        model_hparams["term_use_mpnn"] = False
    if "matches" not in model_hparams.keys():
        model_hparams["matches"] = "resnet"
    if "struct2seq_linear" not in model_hparams.keys():
        model_hparams['struct2seq_linear'] = False
    if "energies_gvp" not in model_hparams.keys():
        model_hparams['energies_gvp'] = False
    if "num_sing_stats" not in model_hparams.keys():
        model_hparams['num_sing_stats'] = 0
    if "num_pair_stats" not in model_hparams.keys():
        model_hparams['num_pair_stats'] = 0
    if "contact_idx" not in model_hparams.keys():
        model_hparams['contact_idx'] = False
    if "fe_dropout" not in model_hparams.keys():
        model_hparams['fe_dropout'] = 0.1
    if "fe_max_len" not in model_hparams.keys():
        model_hparams['fe_max_len'] = 1000
    if "cie_dropout" not in model_hparams.keys():
        model_hparams['cie_dropout'] = 0.1

    if "num_ensembles" in run_hparams.keys():
        model_hparams['num_ensembles'] = run_hparams['num_ensembles']
    else:
        run_hparams['num_ensembles'] = 1
        model_hparams['num_ensembles'] = 1

    if "use_flex" not in model_hparams.keys():
        model_hparams["use_flex"] = False
        model_hparams["flex_type"] = ""

    print(model_hparams)

    # Initialize the model
    terminator = TERMinator(hparams=model_hparams, device=dev)

    # Load weights from the best checkpoint during training
    best_checkpoint_state = torch.load(os.path.join(args.model_dir, 'net_best_checkpoint.pt'), map_location=dev)
    best_checkpoint = best_checkpoint_state['state_dict']
    terminator.load_state_dict(best_checkpoint)
    terminator.to(dev)

    if args.store:
        try:
            os.mkdir('output')
        except OSError as error:
            print(error)

    # Run the model in eval mode, compute the loss and generate energy tables
    scorer = None
    first = True
    print('score structure...')
    start = stop = 0
    for idx,(packaged_structure_data) in enumerate(dataset_iter):
        start = time.time()
        if scorer == None:
            scorer = singlechainScorer(packaged_structure_data,terminator,dev)
        else:
            scorer.set_new_structure(packaged_structure_data)
        scorer.get_scores()
        scorer.write_scores()
        stop = time.time()
        if args.store:
            scorer.pickle_network_output_structure('output/')
        stop = time.time()
        print('Elapsed time:',stop-start)
    
    scorer.close_files()

    print('Done!')
