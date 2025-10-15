import argparse
import json
from multiprocessing.sharedctypes import Value
import os,io,time
import pickle
import numpy as np
import sys
import regex as re
import multiprocessing as mp
import traceback

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

from scripts.design_sequence.ilp_utilities import optimizeEnergyTableViaLP

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

def _raise_error(error):
    """Wrapper for error handling without crashing"""
    traceback.print_exception(Exception, error, None)

def main_mp(args, binder_list):
    print(f"Working with {args.n} processes")
    pool = mp.Pool(args.n)
    batch_size = len(binder_list)//args.n
    binder_list_batches = [binder_list[i:i+batch_size] for i in range(0,len(binder_list),batch_size)]
    for i,binder_list_batch in enumerate(binder_list_batches):
        print(f"Working with subset of {len(binder_list_batch)} structures: {binder_list_batch}")
        res = pool.apply_async(main, args=(args, set(binder_list_batch), i), error_callback=_raise_error)
    pool.close()
    pool.join()

def main_array(args, binder_list):
    batch_size = len(binder_list)//args.n_batches + (len(binder_list) % args.n_batches > 0) # round up
    batch_start = batch_size*args.batch_index
    batch_end = min(batch_size*(args.batch_index+1),len(binder_list))
    print(f"{batch_size} structures split into {args.n_batches} batches")
    print(f"Batch defined as structures: [{batch_start},{batch_end})")
    binder_subset = binder_list[batch_start:batch_end]
    print(f"binder subset has {len(binder_subset)} members: {binder_subset}")
    main(args,binder_subset,args.batch_index)

def get_seq_const_from_name(native_seq,target_chain_len,bridge_name):
    # if this is an extension of a native fragment, constrain that region's sequence 
    # note: this function assumes the region of the native peptide that is to be maintained is at the n-terminus
    native_first_match = re.search("4FXE-ARG81-relax-noHyd-(\d+)-(\d+)__B_(\d)_(\d)_(\d)_",bridge_name)
    native_second_match = re.search("_(\d)_(\d)_(\d)_4FXE-ARG81-relax-noHyd-(\d+)-(\d+)__B",bridge_name)
    native_target_seq = native_seq[:target_chain_len]
    if native_first_match:
        # c-terminal extension
        native_peptide_fragment_len = (int(native_first_match[2]) - int(native_first_match[1]) + 1 - int(native_first_match[3])) # final residue number, first residue number, and amount that was resected from native seed
        native_peptide_fragment_seq = native_seq[target_chain_len:target_chain_len+native_peptide_fragment_len]
        return native_target_seq + native_peptide_fragment_seq + (len(native_seq)-len(native_target_seq)-native_peptide_fragment_len)*'X'

    elif native_second_match:
        # nterminal extension
        native_peptide_fragment_len = (int(native_second_match[5]) - int(native_second_match[4]) + 1 - int(native_second_match[3])) # final residue number, first residue number, and amount that was resected from native seed
        native_peptide_fragment_len = native_peptide_fragment_len - 3 # hacky fix, the wrong sequence was copied at the first three residues
        native_peptide_fragment_seq = native_seq[-native_peptide_fragment_len:]
        return native_target_seq + (len(native_seq)-len(native_target_seq)-native_peptide_fragment_len)*'X' + native_peptide_fragment_seq
    else:
        return None


def get_native_seq(native_seq,start_res_idx,end_res_idx):
    return native_seq[start_res_idx:end_res_idx+1]

def main(args, binder_subset, process_number):
    if (args.binder_dataset):
        if args.dev == 'cpu':
            dataset = BinderScoringIterableDataset(args.binder_dataset,args.target_pdb,500,binder_subset,complex_only=True)
        else:
            dataset = BinderScoringIterableDataset(args.binder_dataset,args.target_pdb,27500,binder_subset,complex_only=True)
        dataset_iter = iter(dataset)
    # elif (args.complex_dataset):
    #     dataset = ComplexScoringDataset(args.complex_dataset)
    #     dataset_iter = iter(dataset)
    else:
        raise ValueError("Must provide either --binder_dataset")

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
    terminator = TERMinator(hparams=model_hparams, device=args.dev)

    # Load weights from the best checkpoint during training
    best_checkpoint_state = torch.load(os.path.join(args.model_dir, 'net_best_checkpoint.pt'), map_location=args.dev)
    best_checkpoint = best_checkpoint_state['state_dict']
    terminator.load_state_dict(best_checkpoint)
    terminator.to(args.dev)
    terminator.eval()
    torch.set_grad_enabled(False)

    # if args.store:
    #     try:
    #         os.mkdir('output')
    #     except OSError as error:
    #         print(error)

    optimizeetab = optimizeEnergyTableViaLP(process_number)

    # Run the model in eval mode, generate energy tables, and perform optimization
    first = True
    nbinders = 0
    if (args.binder_dataset):
        print('Design binders against a single target')
        start = stop = 0
        for packaged_complex_data in dataset_iter:
            nbinders += len(packaged_complex_data['ids'])
            print(f"loaded {len(packaged_complex_data['ids'])} binders, {nbinders} total")
            start = time.time()
            
            # run COORDinator and generate a Potts Model for each structure in the batch
            _to_dev(packaged_complex_data, args.dev)
            max_seq_len = max(packaged_complex_data['seq_lens'].tolist())
            try:
                etab, E_idx, _ = terminator(packaged_complex_data, max_seq_len)
                print('etab: ',etab.shape)
                print('E_idx:',E_idx.shape)
            except Exception as e:
                print(packaged_complex_data['ids'])
                raise e
            net_out_list = []
            for idx in range(etab.shape[0]):
                l, k = etab[idx].shape[0:2]

                # etab and E_idx were padded to fit in the batch, need to remove the extra positions
                l_crop = packaged_complex_data['seq_lens'][idx].item()
                # print('l',l)
                # print('l_crop',l_crop)
                # print(idx,'etab_shape',etab[idx].shape)
                # print(idx,'E_idx_shape',E_idx[idx].shape)
                # print('etab after slice',etab[idx][:l_crop,:min(l_crop,k)].shape)
                # print('E_idx after slice',E_idx[idx][:l_crop,:min(l_crop,k)].shape)

                # print('E_idx before slice',E_idx)
                net_out = {
                    'id': packaged_complex_data['ids'][idx],
                    # 'etab': etab[idx].view(l, k, 20, 20).cpu(),
                    # 'E_idx': E_idx[idx].cpu(),
                    'etab': etab[idx][:l_crop,:min(l_crop,k)].view(l_crop, min(l_crop,k), 20, 20).cpu(),
                    'E_idx': E_idx[idx][:l_crop,:min(l_crop,k)].cpu(),
                    'seq':interfaceScorer.tensorToSeq(packaged_complex_data['seqs'][idx][:l_crop].cpu()),
                    'res_info':packaged_complex_data['res_info'][idx],
                    'chain_lens':packaged_complex_data['chain_lens'][idx]
                }
                print('id',net_out['id'])
                print('etab',net_out['etab'].shape)
                print('E_idx',net_out['E_idx'].shape)
                print('seq',len(net_out['seq']))
                print('res_info',len(net_out['res_info']))
                print('chain_lens',len(net_out['chain_lens']))
                net_out_list.append(net_out)
            
            print(f"Ran the network in eval mode with {len(net_out_list)} structures")
            # Now optimize each energy table using ILP
            for net_out in net_out_list:
                print("Given energy tables and constraints, design sequences")
                # define constraint sequence
                native_seq = net_out['seq']
                print(native_seq)
                bridge_name = net_out['id']
                target_chain_len = net_out['chain_lens'][0]
                print(target_chain_len)
                print(native_seq[:target_chain_len])
                target_constraint_seq = native_seq[:target_chain_len] + len(native_seq[target_chain_len:])*'X'
                print('constraint sequence: ',target_constraint_seq)
                
                constraint_seq = get_seq_const_from_name(native_seq,target_chain_len,bridge_name)
                if constraint_seq == None:
                    print("could not generate constraint sequence from name, skipping")
                    continue
                # # if this is an extension of a native fragment, constrain that region's sequence 
                # # match = re.search("(\d+)-(\d+)_.{1,6}_(\d)_(\d)_(\d)_",bridge_name)
                # native_first_match = re.search("4FXE-ARG81-relax-noHyd-(\d+)-(\d+)__B_(\d)_(\d)_(\d)_",bridge_name)
                # native_second_match = re.search("_(\d)_(\d)_(\d)_4FXE-ARG81-relax-noHyd-(\d+)-(\d+)__B",bridge_name)
                # if native_first_match is not None:
                #     end = net_out['chain_lens'][0] + (int(native_first_match[2]) - int(native_first_match[1]) + 1 - int(native_first_match[3])) # final residue number, first residue number, and amount that was resected from native seed
                #     print(end)
                #     print(native_seq[target_chain_len:end])
                #     print(len(native_seq[end:])*'X')
                #     constraint_seq = native_seq[:target_chain_len] + native_seq[target_chain_len:end] + len(native_seq[end:])*'X'
                # elif native_second_match is not None:
                #     end = net_out['chain_lens'][0] + (int(native_second_match[5]) - int(native_second_match[4]) + 1 - int(native_second_match[3])) # final residue number, first residue number, and amount that was resected from native seed
                #     print(end)
                #     print(native_seq[target_chain_len:end])
                #     print(len(native_seq[end:])*'X')
                #     constraint_seq = native_seq[:target_chain_len] + len(native_seq[target_chain_len:end])*'X' +native_seq[target_chain_len+end:]
                    
                assert len(constraint_seq) == len(native_seq)
                print('constraint sequence: ',constraint_seq)

                # # perform optimization (only constrain with target)
                # optimizeetab.knn_etab_to_sequence(net_out['id'],net_out['etab'], net_out['E_idx'], net_out['res_info'], "target", target_constraint_seq, native_seq)

                # perform optimization (constrain with target + native sequence)
                optimizeetab.knn_etab_to_sequence(net_out['id'],net_out['etab'], net_out['E_idx'], net_out['res_info'], "target_relB", constraint_seq, native_seq)

            stop = time.time()
            print('Elapsed time:',stop-start)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser('')
    parser.add_argument('--binder_dataset', help='Multi-entry PDB file containing a target structure and various binder structures')
    # parser.add_argument('--complex_dataset', help='A file where each line is the path to a PDB file describing a complex. NOTE: peptide must come after the protein structure')
    parser.add_argument('--binder_list', help='A file where each line the name of a binder in the binder_dataset file that should be scored', required=True)
    parser.add_argument('--target_pdb',help='pdb file describing the structure of the target protein', default='')
    # parser.add_argument('--fix_native_seq',help='the name of the region of the peptide that has a fixed sequence',action=)
    parser.add_argument('--model_dir', help='trained model folder', required=True)
    parser.add_argument('--dev', help='device to use', default='cuda:0')
    parser.add_argument('--n',help='the number of cores to use', default=1, type=int)
    parser.add_argument('--n_batches',help='the total number of batches to divide the binders into', default=1, type=int)
    parser.add_argument('--batch_index',help='the index of the batch assigned to this process', default=0, type=int)
    # parser.add_argument('--store', help='if given, will write the output to a JSON file', action=argparse.BooleanOptionalAction)
    args = parser.parse_args()
    print(sys.argv)

    if torch.cuda.device_count() == 0:
        args.dev = "cpu"
    print(f"device: {args.dev}")

    # Initialize the special dataloader for binders
    binder_subset = get_binder_names(args.binder_list)
    print(f"Will attempt to load {len(binder_subset)} total structures from the binder dataset")

    if (args.n > 1):
        main_mp(args,list(binder_subset))
    else:
        main_array(args,list(binder_subset))

    print('Done!')
