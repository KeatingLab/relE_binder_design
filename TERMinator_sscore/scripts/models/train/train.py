"""Train TERMinator model.

Usage:
    .. code-block::

        python train.py \\
            --dataset <dataset_dir> \\
            --model_hparams <model_hparams_file_path> \\
            --run_hparams <run_hparams_file_path> \\
            --run_dir <run_dir> \\
            [--train <train_split_file>] \\
            [--validation <val_split_file>] \\
            [--test <test_split_file>] \\
            [--out_dir <out_dir>] \\
            [--dev <device>] \\
            [--epochs <num_epochs>]
            [--lazy]

    If :code:`--out_dir <out_dir>` is not set, :code:`net.out` will be dumped
    into :code:`<run_dir>`.

    For any of the split files, if the option is not provided, :code:`train.py` will
    look for them within :code:`<dataset_dir>`.

See :code:`python train.py --help` for more info.
"""

import argparse
import datetime
import logging
import copy
import json
import os
os.environ["CUDA_LAUNCH_BLOCKING"] = "1"
import pickle
import sys
import re
from threading import local
print(sys.path)
import optuna
from optuna.trial import TrialState
sys.path.insert(0, '/home/gridsan/sswanson/local_code_mirror/TERMinator_sscore/early-stopping-pytorch')
import torch
import torch.nn as nn
from torch.utils.data import DataLoader
from torch.utils.tensorboard import SummaryWriter
import torch.distributed as dist
from torch.nn.parallel import DistributedDataParallel as DDP
from pytorchtools import EarlyStopping

from terminator.data.data import (TERMLazyDataset, TERMBatchSamplerWrapper, TERMDataset, TERMLazyBatchSamplerWrapper)
from terminator.models.TERMinator import TERMinator
from terminator.utils.model.loop_utils import _log_rank_0, run_epoch
from terminator.utils.model.loss_fn import construct_loss_fn

# for autosummary import purposes
# pylint: disable=wrong-import-order,wrong-import-position
sys.path.insert(0, os.path.dirname(__file__))
from terminator.utils.model.default_hparams import DEFAULT_MODEL_HPARAMS, DEFAULT_TRAIN_HPARAMS
from terminator.utils.model.optim import get_std_opt

# pylint: disable=unspecified-encoding

torch.set_printoptions(threshold=10000)
torch.set_printoptions(linewidth=1000)
torch.set_printoptions(precision=2)
torch.set_printoptions(sci_mode=False)
torch.set_printoptions(profile="full")


def _setup_hparams(args, trial, device_rank):
    """ Setup the hparams dictionary using defaults and return it

    Args
    ----
    args : argparse.Namespace
        Parsed arguments
    trial : optuna.Trial
        Optuna process for hyperparameter optimization
    

    Returns
    -------
    model_hparams : dict
        Fully configured model hparams dictionary (see :code:`terminator/utils/model/default_hparams.py`)
    run_hparams : dict
        Fully configured training run hparams dictionary (see :code:`terminator/utils/model/default_hparams.py`)
    """
    def _load_hparams(hparam_path, default_hparams, output_name):
        print("loading params")
        # load hparams
        hparams_path = os.path.join(args.run_dir, output_name)
        if device_rank > 0 and os.path.isfile(hparams_path) and (not os.path.getsize(hparams_path) == 0) and (args.n_trials == 0):
            hparams = json.load(open(hparams_path, 'r'))
            for hparam in default_hparams:
                if hparam not in hparams:
                    hparams[hparam] = default_hparams[hparam]
            return hparams, None

        hparams = json.load(open(hparam_path, 'r'))
        optuna_hparams = {}
        for key, default_val in default_hparams.items():
            if key not in hparams:
                hparams[key] = default_val
        hparams_copy = copy.deepcopy(hparams)
        for key, default_val in hparams_copy.items():
            if re.search('^optuna_', key) != None:
                key_split = key.split('_')
                key = "_".join(key_split[2:])
                print(key)
                optuna_type = key_split[1]
                optuna_vals = default_val.split(":")
                if len(optuna_vals) > 2:
                    step = optuna_vals[1]
                    optuna_vals[1] = optuna_vals[2]
                else:
                    step = None
                if optuna_type == 'categorical':
                    hparams[key] = trial.suggest_categorical(key, optuna_vals)
                elif optuna_type == 'loguniform':
                    hparams[key] = trial.suggest_loguniform(key, float(optuna_vals[0]), float(optuna_vals[1]))
                elif optuna_type == 'float':
                    print(key, optuna_vals[0], optuna_vals[1], step)
                    hparams[key] = trial.suggest_float(key, float(optuna_vals[0]), float(optuna_vals[1]), step=float(step))
                    hparams[key] = 0.02 ## JFM
                else:
                    hparams[key] = trial.suggest_int(key, int(optuna_vals[0]), int(optuna_vals[1]), step=int(step))
                optuna_hparams[key] = hparams[key]

        return hparams, optuna_hparams
    
    def _save_hparams(hparams, output_name, check_conflict = True):
        hparams_path = os.path.join(args.run_dir, output_name)
        if check_conflict and os.path.isfile(hparams_path) and (not os.path.getsize(hparams_path) == 0):
            previous_hparams = json.load(open(hparams_path, 'r'))
            if previous_hparams != hparams:
                raise Exception('Given hyperparameters do not agree with previous hyperparameters.')
        else:
            json.dump(hparams, open(hparams_path, 'w'))


    model_hparams, optuna_model_hparams = _load_hparams(args.model_hparams, DEFAULT_MODEL_HPARAMS, 'model_hparams.json')
    run_hparams, optuna_run_hparams = _load_hparams(args.run_hparams, DEFAULT_TRAIN_HPARAMS, 'run_hparams.json')
    if args.base_run_dir:
        args.run_dir = args.base_run_dir
    if optuna_model_hparams is not None and optuna_run_hparams is not None:
        optuna_hparams = dict(optuna_model_hparams, **optuna_run_hparams)
        args.base_run_dir = args.run_dir
        args.run_dir = args.run_dir if not optuna_hparams else args.run_dir + "/run_" + "_".join([param + "_" + str(optuna_hparams[param]) for param in optuna_hparams])
    if not os.path.isdir(args.run_dir):
        os.makedirs(args.run_dir, exist_ok=True)
    if (not dist.is_initialized()) or (dist.get_rank() == 0):
        _save_hparams(model_hparams, 'model_hparams.json')
        _save_hparams(run_hparams, 'run_hparams.json')

    return model_hparams, run_hparams


def _setup_dataloaders(args, run_hparams, model_hparams):
    """ Setup dataloaders needed for training

    Args
    ----
    args : argparse.Namespace
        Parsed arguments
    run_hparams : dict
        Fully configured hparams dictionary (see :code:`terminator/utils/model/default_hparams.py`)

    Returns
    -------
    train_dataloader, val_dataloader, test_dataloader : torch.utils.data.DataLoader
        DataLoaders for the train, validation, and test datasets
    """
    kwargs = {}
    kwargs['num_workers'] = 16

    # set up dataloaders
    train_ids = []
    with open(args.train, 'r') as f:
        for line in f:
            train_ids += [line.strip()]
    validation_ids = []
    with open(args.validation, 'r') as f:
        for line in f:
            validation_ids += [line.strip()]
    test_ids = []
    with open(args.test, 'r') as f:
        for line in f:
            test_ids += [line.strip()]
    if args.lazy:
        train_dataset = TERMLazyDataset(args.dataset, pdb_ids=train_ids, min_protein_len=model_hparams['min_protein_len'])
        val_dataset = TERMLazyDataset(args.dataset, pdb_ids=validation_ids, min_protein_len=model_hparams['min_protein_len'])
        test_dataset = TERMLazyDataset(args.dataset, pdb_ids=test_ids, min_protein_len=model_hparams['min_protein_len'])
        train_batch_sampler = TERMLazyBatchSamplerWrapper(ddp=dist.is_initialized())
        train_batch_sampler = train_batch_sampler.sampler(train_batch_sampler.ddp, train_dataset,
                                                   batch_size=run_hparams['train_batch_size'],
                                                   shuffle=run_hparams['shuffle'],
                                                   semi_shuffle=run_hparams['semi_shuffle'],
                                                   sort_data=run_hparams['sort_data'],
                                                   term_matches_cutoff=run_hparams['term_matches_cutoff'],
                                                   max_term_res=run_hparams['max_term_res'],
                                                   max_seq_tokens=run_hparams['max_seq_tokens'],
                                                   term_dropout=run_hparams['term_dropout'],
                                                   noise_type=run_hparams['noise_type'],
                                                   noise=run_hparams['noise'])
        if 'test_term_matches_cutoff' in run_hparams:
            test_term_matches_cutoff = run_hparams['test_term_matches_cutoff']
        else:
            test_term_matches_cutoff = run_hparams['term_matches_cutoff']
        val_batch_sampler = TERMLazyBatchSamplerWrapper(ddp=dist.is_initialized())
        val_batch_sampler = val_batch_sampler.sampler(val_batch_sampler.ddp, val_dataset,
                                                 batch_size=1,
                                                 shuffle=False,
                                                 term_matches_cutoff=test_term_matches_cutoff)
        test_batch_sampler = TERMLazyBatchSamplerWrapper(ddp=dist.is_initialized())
        test_batch_sampler = test_batch_sampler.sampler(test_batch_sampler.ddp, test_dataset,
                                                  batch_size=1,
                                                  shuffle=False,
                                                  term_matches_cutoff=test_term_matches_cutoff)
    else:
        train_dataset = TERMDataset(args.dataset, run_hparams['pdb_dataset'], run_hparams['flex_folder'],
                                    model_hparams['flex_type'], run_hparams['noise_level'], run_hparams['bond_length_noise_level'], pdb_ids=train_ids,
                                    max_protein_len=run_hparams['max_seq_tokens'], num_ensembles=run_hparams['num_ensembles'], min_protein_len=model_hparams['min_protein_len'])
        val_dataset = TERMDataset(args.dataset, run_hparams['pdb_dataset'], run_hparams['flex_folder'],
                                    model_hparams['flex_type'], pdb_ids=validation_ids,
                                    max_protein_len=run_hparams['max_seq_tokens'], num_ensembles=run_hparams['num_ensembles'], min_protein_len=model_hparams['min_protein_len'])
        test_dataset = TERMDataset(args.dataset, run_hparams['pdb_dataset'], run_hparams['flex_folder'],
                                    model_hparams['flex_type'], pdb_ids=test_ids,
                                    max_protein_len=run_hparams['max_seq_tokens'], num_ensembles=run_hparams['num_ensembles'], min_protein_len=model_hparams['min_protein_len'])

        train_batch_sampler = TERMBatchSamplerWrapper(ddp=dist.is_initialized())
        train_batch_sampler = train_batch_sampler.sampler(train_batch_sampler.ddp, train_dataset, args.dev,
                                               batch_size=run_hparams['train_batch_size'],
                                               shuffle=run_hparams['shuffle'],
                                               semi_shuffle=run_hparams['semi_shuffle'],
                                               sort_data=run_hparams['sort_data'],
                                               max_term_res=run_hparams['max_term_res'],
                                               max_seq_tokens=run_hparams['max_seq_tokens'],
                                               flex_type=model_hparams['flex_type'],
                                               noise_level=run_hparams['noise_level'],
                                               bond_length_noise_level=run_hparams['bond_length_noise_level'],
                                               num_ensembles=run_hparams['num_ensembles'])
        val_batch_sampler = TERMBatchSamplerWrapper(ddp=dist.is_initialized())
        val_batch_sampler = val_batch_sampler.sampler(val_batch_sampler.ddp, val_dataset, args.dev, batch_size=1, shuffle=False, flex_type=model_hparams['flex_type'], num_ensembles=run_hparams['num_ensembles'])
        test_batch_sampler = TERMBatchSamplerWrapper(ddp=dist.is_initialized())
        test_batch_sampler = test_batch_sampler.sampler(test_batch_sampler.ddp, test_dataset, args.dev, batch_size=1, shuffle=False, flex_type=model_hparams['flex_type'], num_ensembles=run_hparams['num_ensembles'])

    train_dataloader = DataLoader(train_dataset,
                                  batch_sampler=train_batch_sampler,
                                  collate_fn=train_batch_sampler.package,
                                  pin_memory=True,
                                  **kwargs)
    val_dataloader = DataLoader(val_dataset,
                                batch_sampler=val_batch_sampler,
                                collate_fn=val_batch_sampler.package,
                                pin_memory=True,
                                **kwargs)
    test_dataloader = DataLoader(test_dataset,
                                 batch_sampler=test_batch_sampler,
                                 collate_fn=test_batch_sampler.package,
                                 **kwargs)

    return train_dataloader, val_dataloader, test_dataloader, train_batch_sampler


def _load_checkpoint(run_dir, dev, finetune=False):
    """ If a training checkpoint exists, load the checkpoint. Otherwise, setup checkpointing initial values.

    Args
    ----
    run_dir : str
        Path to directory containing the training run checkpoint, as well the tensorboard output.

    Returns
    -------
    dict
        Dictionary containing
        - "best_checkpoint_state": the best checkpoint state during the run
        - "last_checkpoint_state": the most recent checkpoint state during the run
        - "best_checkpoint": the best model parameter set during the run
        - "best_validation": the best validation loss during the run
        - "last_optim_state": the most recent state of the optimizer
        - "start_epoch": what epoch to resume training from
        - "writer": SummaryWriter for tensorboard
        - "training_curves": pairs of (train_loss, val_loss) representing the training and validation curves
    """

    if os.path.isfile(os.path.join(run_dir, 'net_best_checkpoint.pt')):
        best_checkpoint_state = torch.load(os.path.join(run_dir, 'net_best_checkpoint.pt'), map_location=torch.device(dev))
        last_checkpoint_state = torch.load(os.path.join(run_dir, 'net_last_checkpoint.pt'), map_location=torch.device(dev))
        best_checkpoint = best_checkpoint_state['state_dict']
        best_validation = best_checkpoint_state['val_loss']
        last_optim_state = last_checkpoint_state["optimizer_state"]
        start_epoch = last_checkpoint_state['epoch'] + 1
        writer = SummaryWriter(log_dir=os.path.join(run_dir, 'tensorboard'), purge_step=start_epoch + 1)
        training_curves = last_checkpoint_state["training_curves"]
    else:
        best_checkpoint_state, last_checkpoint_state = None, None
        best_checkpoint = None
        best_validation = 10e8
        last_optim_state = None
        start_epoch = 0
        writer = SummaryWriter(log_dir=os.path.join(run_dir, 'tensorboard'))
        training_curves = {"train_loss": [], "val_loss": []}
        if finetune: # load existing model for finetuning
            best_checkpoint_state = torch.load(os.path.join(run_dir, 'net_original.pt'), map_location=torch.device(dev))
            best_checkpoint = best_checkpoint_state['state_dict']

    return {"best_checkpoint_state": best_checkpoint_state,
            "last_checkpoint_state": last_checkpoint_state,
            "best_checkpoint": best_checkpoint,
            "best_validation": best_validation,
            "last_optim_state": last_optim_state,
            "start_epoch": start_epoch,
            "writer": writer,
            "training_curves": training_curves}


def _setup_model(model_hparams, run_hparams, checkpoint, dev, local_rank):
    """ Setup a TERMinator model using hparams, a checkpoint if provided, and a computation device.

    Args
    ----
    model_hparams : dict
        Fully configured model hparams dictionary (see :code:`terminator/utils/model/default_hparams.py`)
    run_hparams : dict
        Fully configured training run hparams dictionary (see :code:`terminator/utils/model/default_hparams.py`)
    checkpoint : OrderedDict or None
        Model parameters
    dev : str
        Computation device to use
    local_rank : int
        Indicator of local process rank for ddp

    Returns
    -------
    terminator : TERMinator or nn.DataParallel(TERMinator)
        Potentially parallelized TERMinator to use for training
    terminator_module : TERMinator
        Inner TERMinator, unparallelized
    """
    model_hparams['num_ensembles'] = run_hparams['num_ensembles']
    terminator = TERMinator(hparams=model_hparams, device=dev)
    if checkpoint is not None:
        missing_keys, unexpected_keys = terminator.load_state_dict(checkpoint,strict=False)
        if len(unexpected_keys) > 0:
            print("Checkpoint has unexpected keys:",unexpected_keys)
            # raise ValueError("Checkpoint has unexpected keys:",unexpected_keys)
        finetune = ('finetune_nlcpl' in run_hparams and run_hparams['finetune_nlcpl'] == True) or ('finetune_sscore' in run_hparams and run_hparams['finetune_sscore'] == True)
        if not finetune and len(missing_keys) > 0:
            raise ValueError("Checkpoint has missing keys and not set to finetune:",missing_keys)

    _log_rank_0(terminator)
    _log_rank_0("terminator hparams: " + str(terminator.hparams))

    if torch.cuda.device_count() > 1 and dev != "cpu" and not dist.is_initialized():
        terminator = nn.DataParallel(terminator)
        terminator_module = terminator.module
    else:
        terminator_module = terminator
    terminator.to(dev)
    if dist.is_initialized():
        terminator = DDP(
                terminator,
                device_ids=[local_rank],
                output_device=local_rank,
                find_unused_parameters=True
        )
        _log_rank_0("DDP setup finished")

    return terminator, terminator_module


def out_results(training_curves, terminator, terminator_module, run_dir, best_checkpoint, test_dataloader, loss_fn, dev, args, device_rank):
    """ Evaluate best checkpoint on test set and output results.

    Args
    ----
    training_curves : dict
        Training set and validation set losses for all epochs
    terminator : TERMinator or nn.DataParallel(TERMinator)
        Potentially parallelized TERMinator to use for training
    run_dir : str
        Path to directory containing the training run checkpoint, as well the tensorboard output.
    best_checkpoint : OrderedDict or None
        Best model parameters
    test_dataloader : torch.utils.data.DataLoader
        DataLoaders for the train, validation, and test datasets.
    loss_fn : function
        The constructed loss function.
    dev : str
        Computation device to use
    args : argparse.Namespace
        Parsed arguments
    Returns
    -------
    no return value
    """
    # save model params
    _log_rank_0(str(training_curves))
    if device_rank < 1:
        torch.save(terminator_module.state_dict(), os.path.join(run_dir, 'net_last.pt'))
        torch.save(best_checkpoint, os.path.join(run_dir, 'net_best.pt'))

    # test
    terminator_module.load_state_dict(best_checkpoint)
    test_loss, test_ld, dump = run_epoch(terminator, test_dataloader, loss_fn, grad=False, test=True, dev=dev)
    _log_rank_0(f"test loss {test_loss} test loss dict {test_ld}")
    # dump outputs
    if args.out_dir and device_rank < 1:
        if not os.path.isdir(args.out_dir):
            os.mkdir(args.out_dir)
        net_out_path = os.path.join(args.out_dir, "net.out")
    else:
        net_out_path = os.path.join(run_dir, "net.out")
    # save etab outputs for dTERMen runs
    with open(net_out_path, 'wb') as fp:
        pickle.dump(dump, fp)
    return test_loss


def objective_func(args, run_hparams, model_hparams, local_rank, device_rank, trial):
    """ Train TERMinator """
    dev = args.dev
    
    # setup dataloaders
    train_dataloader, val_dataloader, test_dataloader, train_batch_sampler = _setup_dataloaders(args, run_hparams, model_hparams)
    # load checkpoint
    checkpoint_dict = _load_checkpoint(args.run_dir, dev, run_hparams['finetune_nlcpl'] or run_hparams['finetune_sscore'])
    best_validation = checkpoint_dict["best_validation"]
    best_checkpoint = checkpoint_dict["best_checkpoint"]
    start_epoch = checkpoint_dict["start_epoch"]
    last_optim_state = checkpoint_dict["last_optim_state"]
    writer = checkpoint_dict["writer"]
    training_curves = checkpoint_dict["training_curves"]

    isDataParallel = True if torch.cuda.device_count() > 1 and dev != "cpu" else False
    finetune_nlcpl = run_hparams["finetune_nlcpl"]
    finetune_sscore = run_hparams["finetune_sscore"]

    # construct terminator, loss fn, and optimizer
    terminator, terminator_module = _setup_model(model_hparams, run_hparams, best_checkpoint, dev, local_rank)
    loss_fn = construct_loss_fn(run_hparams)
    optimizer = get_std_opt(terminator.parameters(),
                            d_model=model_hparams['energies_hidden_dim'],
                            regularization=run_hparams['regularization'],
                            state=last_optim_state,
                            finetune=finetune_nlcpl or finetune_sscore,
                            finetune_lr=run_hparams["finetune_lr"])

    # initialize early stopping
    early_stopping = EarlyStopping(patience=run_hparams['patience'], verbose=True)
    try:
        _log_rank_0('iterating')
        for epoch in range(start_epoch, args.epochs):
            _log_rank_0('epoch ' + str(epoch))
            if dist.is_initialized():
                train_batch_sampler.set_epoch(epoch)

            epoch_loss, epoch_ld, _ = run_epoch(terminator, train_dataloader, loss_fn, optimizer=optimizer, grad=True, dev=dev, finetune_nlcpl=finetune_nlcpl, finetune_sscore=finetune_sscore, isDataParallel=isDataParallel)
            _log_rank_0('epoch loss: ' + str(epoch_loss) + ', epoch_ld: ' + str(epoch_ld))
            if device_rank < 1:
                writer.add_scalar('training loss', epoch_loss, epoch)

            # validate
            val_loss, val_ld, _ = run_epoch(terminator, val_dataloader, loss_fn, grad=False, dev=dev)
            _log_rank_0('val loss: ' + str(val_loss) + ', val ld: ' + str(val_ld))
            writer.add_scalar('val loss', val_loss, epoch)

            training_curves["train_loss"].append((epoch_loss, epoch_ld))
            training_curves["val_loss"].append((val_loss, val_ld))
            # comp = (val_ld['sortcery_loss']['loss'] < best_validation) if finetune else (val_loss < best_validation)

            # save a state checkpoint
            checkpoint_state = {
                'epoch': epoch,
                'state_dict': terminator_module.state_dict(),
                'best_model': (val_loss < best_validation),
                'val_loss': best_validation,
                'optimizer_state': optimizer.state_dict(),
                'training_curves': training_curves
            }
            if device_rank < 1:
                torch.save(checkpoint_state, os.path.join(args.run_dir, 'net_last_checkpoint.pt'))
            if (val_loss < best_validation): # if comp:
                # if finetune:
                #     best_validation = val_ld['sortcery_loss']['loss']
                # else:
                best_validation = val_loss
                best_checkpoint = copy.deepcopy(terminator_module.state_dict())
                if device_rank < 1:
                    torch.save(checkpoint_state, os.path.join(args.run_dir, 'net_best_checkpoint.pt'))

            # Check for early stopping
            early_stopping(val_loss, terminator)
    
            if early_stopping.early_stop:
                _log_rank_0(f"Early stopping at epoch {epoch}")
                test_loss = out_results(training_curves, terminator, terminator_module, args.run_dir, best_checkpoint, test_dataloader, loss_fn, dev, args, device_rank)
                if args.n_trials > 0:
                    raise optuna.exceptions.TrialPruned()
                return test_loss
            if args.n_trials > 0:
                trial.report(val_loss, step=epoch)
                if trial.should_prune():
                    raise optuna.exceptions.TrialPruned()
            if dist.is_initialized() and local_rank != -1:
                dist.barrier()

    except KeyboardInterrupt:
        pass

    test_loss = out_results(training_curves, terminator, terminator_module, args.run_dir, best_checkpoint, test_dataloader, loss_fn, dev, args, device_rank)

    writer.close()
    return test_loss


def main(base_trial):
    if device_rank > 0:
        trial = None
    if args.n_trials > 0 and local_rank > -1:
        trial = optuna.integration.TorchDistributedTrial(base_trial, local_rank) 
    else:
        trial = None
    model_hparams, run_hparams = _setup_hparams(args, trial, device_rank)
    if run_hparams['flex_folder']:
        args.train = os.path.join(run_hparams['flex_folder'], 'train.in')
        args.validation = os.path.join(run_hparams['flex_folder'], 'validation.in')
        args.test = os.path.join(run_hparams['flex_folder'], 'test.in')
    # Check duplication and skip if it's detected.
    if device_rank == 0:
        for t in base_trial.study.trials:
            if t.state != optuna.trial.TrialState.COMPLETE:
                continue

            if t.params == base_trial.params:
                return t.value  # Return the previous value without re-evaluating it.
    test_loss = objective_func(args, run_hparams, model_hparams, local_rank, device_rank, trial)
    return test_loss


if __name__ == '__main__':

    # initialize COORDinator args
    parser = argparse.ArgumentParser('Train COORDinator!')
    parser.add_argument('--dataset', help='input folder .features files in proper directory structure.', required=True)
    parser.add_argument('--model_hparams', help='file path for model hparams', required=True)
    parser.add_argument('--run_hparams', help='file path for run hparams', required=True)
    parser.add_argument('--run_dir', help='path to place folder to store model files', required=True)
    parser.add_argument('--base_run_dir', help='base path to store model files when running with optuna', default=None, type=str)
    parser.add_argument('--train', help='file with training dataset split')
    parser.add_argument('--validation', help='file with validation dataset split')
    parser.add_argument('--test', help='file with test dataset split')
    parser.add_argument('--out_dir',
                        help='path to place test set eval results (e.g. net.out). If not set, default to --run_dir')
    parser.add_argument('--dev', help='device to train on', default='cuda:0')
    parser.add_argument('--epochs', help='number of epochs to train for', default=100, type=int)
    parser.add_argument('--lazy', help="use lazy data loading", action='store_true')
    parser.add_argument('--n_nodes', help="number of cores for use in ddp", default=1, type=int)
    parser.add_argument('--n_trials', help="number of trials for optuna optimization", default=0, type=int)
    parser.add_argument('--optuna_verbose', help="save best model checkpoints for all trials", default=False, type=bool)
    parser.add_argument("--backend", help="Backend for DDP", type=str, default="gloo")
    parsed_args = parser.parse_args()
    args = parsed_args
    # by default, if no splits are provided, read the splits from the dataset folder
    if parsed_args.train is None:
        parsed_args.train = os.path.join(parsed_args.dataset, 'train.in')
    if parsed_args.validation is None:
        parsed_args.validation = os.path.join(parsed_args.dataset, 'validation.in')
    if parsed_args.test is None:
        parsed_args.test = os.path.join(parsed_args.dataset, 'test.in')
    print(parsed_args.train, parsed_args.validation, parsed_args.test)

    # setup ddp
    if args.dev == "cpu":
        local_rank = -1
    else:
        local_rank = int(os.environ["LOCAL_RANK"])
    if parsed_args.n_nodes > 1 and local_rank != -1:
        dist.init_process_group(backend=parsed_args.backend, init_method='env://', timeout=datetime.timedelta(hours=3))
        torch.cuda.set_device(local_rank)
        torch.backends.cudnn.benchmark = True
    if dist.is_initialized():
        args.dev = torch.device("cuda")
        print(f"Local rank: {local_rank}")
        print(f"Device rank: {dist.get_rank()}")
        print(f"World size: {dist.get_world_size()}")
        device_rank = dist.get_rank()
    else:
        device_rank = -1
    # initialize optuna 
    if parsed_args.n_trials > 0:
        study = None
        study = optuna.create_study(direction="minimize", pruner=optuna.pruners.NopPruner())
        if device_rank == 0:
            
            study.optimize(main, n_trials = parsed_args.n_trials)
        else:
            i = -1
            while parsed_args.n_trials > len(set(str(t.params) for t in study.trials)):
                i += 1
                try:
                    print(f"starting trial {i}")
                    main(None)
                except optuna.TrialPruned:
                    pass
            
        
        if device_rank == 0:
            assert study is not None
            pruned_trials = study.get_trials(deepcopy=False, states=[TrialState.PRUNED])
            complete_trials = study.get_trials(deepcopy=False, states=[TrialState.COMPLETE])

            print("Study statistics: ")
            print("  Number of finished trials: ", len(study.trials))
            print("  Number of pruned trials: ", len(pruned_trials))
            print("  Number of complete trials: ", len(complete_trials))

            print("Best trial:")
            trial = study.best_trial

            print("  Value: ", trial.value)

            print("  Params: ")
            for key, value in trial.params.items():
                print("    {}: {}".format(key, value))    
    else:
        main(None)
