"""Score de-novo designed protein-binding peptide backbones using TERMinator energy tables

# The resulting evaluated proteins will be dumped in :code:`<output_dir>` via
# a pickle file :code:`net.out`.

Usage:
    .. code-block::

        python scoreBinders.py \\
            # --dataset <dataset_dir> \\
            # --model_dir <trained_model_dir> \\
            # --output_dir <output_dir> \\
            # [--subset <data_subset_file>] \\
            # [--dev <device>]

    If :code:`subset` is not provided, the entire dataset :code:`dataset` will
    be evaluated.

See :code:`python scoreBinders.py --help` for more info.
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
from terminator.data.data import BinderScoringIterableDataset,ComplexScoringDataset
from terminator.models.TERMinator import TERMinator
from terminator.utils.model.loop_utils import run_epoch,_to_dev
from terminator.utils.model.loss_fn import construct_loss_fn
from terminator.utils.common import int_to_3lt_AA,AA_to_int,int_to_AA
from terminator.utils.model.default_hparams import DEFAULT_MODEL_HPARAMS, DEFAULT_TRAIN_HPARAMS

# from python.score_binders.hbond_utils import bbHBond
from python.score_binders.score_utils import *

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
    parser = argparse.ArgumentParser('Score peptide-binder backbones using TERMinator energies')
    parser.add_argument('--binder_dataset', help='Multi-entry PDB file containing a target structure and various binder structures')
    parser.add_argument('--complex_dataset', help='A file where each line is the path to a PDB file describing a complex. NOTE: peptide must come after the protein structure')
    parser.add_argument('--binder_list', help='A file where each line the name of a binder in the binder_dataset file that should be scored', default='')
    parser.add_argument('--target_pdb',help='pdb file describing the structure of the target protein', default='')
    parser.add_argument('--model_dir', help='trained model folder', required=True)
    parser.add_argument('--seq_mode', help='controls which sequence is used when scoring the binder', default='consensus_aa')
    parser.add_argument('--score_mode', help='controls which pair energies are used when scoring', default='')
    parser.add_argument('--custom_pep_reference', help='', default='')
    parser.add_argument('--target_chain_id', help='In complex mode, optionally pass the protein chain ID(s) (e.g.: "ABC")',default='',type=str)
    parser.add_argument('--binder_chain_id', help='In complex mode, optionally pass the peptide chain ID',default='')
    parser.add_argument('--dev', help='device to use', default='cuda:0')
    parser.add_argument('--store', help='if given, will write the output to a JSON file', action=argparse.BooleanOptionalAction)
    parser.add_argument('--hbond', help='if given, will print out hydrogen bonding information', action=argparse.BooleanOptionalAction)
    args = parser.parse_args()
    print(sys.argv)

    '''
    If scoring many structures against a single target, use `BinderScoringIterableDataset` class. This
    loads from a single multi-entry PDB file, where the first section is the target and the remaining 
    sections are binder structures. A file of this format is written by "generateSeeds" in the interface Generator repo

    If scoring complexes, provide a path to a file where each line is the path to a complex structure.
    The name of the file must be of the format NAME_X_Y, where Y is the single character binder chain ID.
    '''

    if args.hbond:
        raise ValueError('hydrogen bond detection not yet updated to with batched dataloaders')

    dev = args.dev
    if torch.cuda.device_count() == 0:
        dev = "cpu"

    # Initialize the special dataloader for binders
    binder_subset = get_binder_names(args.binder_list)
    n_binders = len(binder_subset) if (binder_subset is not None) else "all"

    # If a simple peptide reference state is being used, no need to load binder structures separately
    load_mode = "complex_only" if args.custom_pep_reference != '' else "binder_and_complex"
    # load_mode = 'binder_and_complex'
    print(f"Will attempt to load {n_binders} total structures from the binder dataset in {load_mode} mode")
    if (args.binder_dataset):
        # dataset = BinderScoringIterableDataset(args.binder_dataset,args.target_pdb)
        dataset = BinderScoringIterableDataset(filename=args.binder_dataset,target_pdb_path=args.target_pdb,max_res_num=27500,binder_subset=binder_subset,mode=load_mode,skip_package=False)
        dataset_iter = iter(dataset)
    elif (args.complex_dataset):
        dataset = ComplexScoringDataset(args.complex_dataset,targetChainID=args.target_chain_id,binderChainID=args.binder_chain_id,mode=load_mode)
        dataset_iter = iter(dataset)
    else:
        raise ValueError("Must provide either --binder_dataset or --complex_dataset")
    # dataloader = DataLoader(dataset,batch_size=1)

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
    binderSubsetNames = get_binder_names(args.binder_list)

    custom_pep_ref_aa_probs = None
    if args.custom_pep_reference != "":
        with open(args.custom_pep_reference,"rb") as file:
            custom_pep_ref_aa_probs = pickle.load(file)
        print(f"Loaded custom peptide reference aa probabilites: {custom_pep_ref_aa_probs}")

    if args.store:
        try:
            os.mkdir('output')
        except OSError as error:
            print(error)

    # Run the model in eval mode, compute the loss and generate energy tables
    scorer = None
    first = True
    if (args.binder_dataset):
        print('score binders against a single target: ',args.score_mode)
        start = stop = 0
        for packaged_binder_complex_data in dataset_iter:
            start = time.time()
            if scorer == None:
                scorer = interfaceScorer(dataset.target_packaged_data,packaged_binder_complex_data,terminator,custom_pep_ref_aa_probs,dev)
                # if args.hbond:
                #     hbDetector = bbHBond(dataset.target_data['coords'])
                #     hbDetector.setBinderResidues(dataset.current_binder_data['coords'])
            else:
                scorer.set_new_binder_and_complex(packaged_binder_complex_data)
                # if args.hbond:
                #     hbDetector.setBinderResidues(dataset.current_binder_data['coords'])
            scorer.get_binder_scores(args.seq_mode,args.score_mode)
            scorer.write_structure_scores(args.score_mode)
            scorer.write_residue_scores(args.score_mode)
            scorer.write_pep_seqs()
            # scorer.write_pep_seq_probs()
            stop = time.time()
            # if args.hbond:
            #     hbDetector.findHBonds()
            #     hbDetector.writeHBonds(scorer.complex_netout['id'],scorer.prot_res_info,scorer.pep_res_info)
            #     hbDetector.writeExposed(scorer.complex_netout['id'],scorer.pep_res_info)
            #     hbDetector.writeHBondsTotalE(scorer.complex_netout['id'])
            #     # hbDetector.reportHBonds(scorer.complex_netout['id'])
            if args.store:
                if first:
                    # only need to store once
                    scorer.pickle_network_output_target('output/')
                    first = False
                scorer.pickle_network_output_binder('output/')
                scorer.pickle_network_output_complex('output/')
            print('Elapsed time:',stop-start)
    else:
        # Score complexes
        print('score complexes: ',args.score_mode)
        start = stop = 0
        for idx,(packaged_target_data,packaged_binder_complex_data) in enumerate(dataset_iter):
            start = time.time()
            if scorer == None:
                scorer = interfaceScorer(packaged_target_data,packaged_binder_complex_data,terminator,custom_pep_ref_aa_probs,dev)
                # if args.hbond:
                #     hbDetector = bbHBond(dataset.target_data['coords'])
            else:
                scorer.set_new_target(packaged_target_data)
                scorer.set_new_binder_and_complex(packaged_binder_complex_data)
                # if args.hbond:
                #     hbDetector.setTargetResidues(dataset.target_data['coords'])
            if binderSubsetNames is not None and scorer.get_pdb_name() not in binderSubsetNames:
                continue
            # else:
            #     print(f"scoring binder {scorer.get_pdb_name(idx)}")
            scorer.get_binder_scores(args.seq_mode,args.score_mode)
            scorer.write_structure_scores(args.score_mode)
            scorer.write_residue_scores(args.score_mode)
            scorer.write_pep_seqs()
            scorer.write_pep_seq_probs()
            stop = time.time()
            # if args.hbond:
            #     hbDetector.setBinderResidues(dataset.complex_data['coords'][scorer.pep_res_range[0]:scorer.pep_res_range[1]])
            #     hbDetector.findHBonds()
            #     hbDetector.writeHBonds(scorer.complex_netout['id'],scorer.prot_res_info,scorer.pep_res_info)
            #     hbDetector.writeExposed(scorer.complex_netout['id'],scorer.pep_res_info)
            #     hbDetector.writeHBondsTotalE(scorer.complex_netout['id'])
            #     # hbDetector.reportHBonds(scorer.complex_netout['id'])
            if args.store:
                scorer.pickle_network_output_target('output/')
                scorer.pickle_network_output_binder('output/')
                scorer.pickle_network_output_complex('output/')
            stop = time.time()
            print('Elapsed time:',stop-start)
    
    scorer.close_files()
    # if args.hbond:
    #     hbDetector.close_files()

    print('Done!')
