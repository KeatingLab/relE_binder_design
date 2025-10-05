"""
Modification version of https://github.com/optuna/optuna/pull/2303 with nccl backend
Optuna example that optimizes multi-layer perceptrons using PyTorch distributed.
In this example, we optimize the validation accuracy of hand-written digit recognition using
PyTorch distributed data parallel and MNIST. We optimize the neural network architecture as well
as the optimizer configuration. As it is too time consuming to use the whole MNIST dataset, we
here use a small subset of it.
You can execute this example with mpirun command as follows:
    $  python -m torch.distributed.launch --nproc_per_node=2 pytorch_distributed_simple.py
Please note that you need to install PyTorch from source if you switch the communication backend
of torch.distributed to "mpi". Please refer to the following document for further details:
https://pytorch.org/tutorials/intermediate/dist_tuto.html#communication-backends
"""

import argparse
import os

import torch
import torch.distributed as dist
import torch.nn as nn
import torch.nn.functional as F
from torch.nn.parallel import DistributedDataParallel as DDP
import torch.optim as optim
import torch.utils.data
from torchvision import datasets
from torchvision import transforms

import optuna


BATCHSIZE = 128
CLASSES = 10
DIR = os.getcwd()
EPOCHS = 1
LOG_INTERVAL = 10
N_TRAIN_EXAMPLES = BATCHSIZE * 30
N_VALID_EXAMPLES = BATCHSIZE * 10


def define_model(trial):
    # We optimize the number of layers, hidden units and dropout ratio in each layer.
    n_layers = trial.suggest_int("n_layers", 1, 3)
    layers = []

    in_features = 28 * 28
    for i in range(n_layers):
        print("about to suggest int")
        out_features = trial.suggest_int("n_units_l{}".format(i), 4, 128)
        print("int suggested")
        layers.append(nn.Linear(in_features, out_features))
        layers.append(nn.ReLU())
        p = trial.suggest_float("dropout_l{}".format(i), 0.2, 0.5)
        layers.append(nn.Dropout(p))

        in_features = out_features
    layers.append(nn.Linear(in_features, CLASSES))
    layers.append(nn.LogSoftmax(dim=1))

    return nn.Sequential(*layers)


def get_mnist():
    # Load MNIST dataset.
    train_dataset = datasets.MNIST(DIR, train=True, transform=transforms.ToTensor())
    train_dataset = torch.utils.data.Subset(train_dataset, indices=range(N_TRAIN_EXAMPLES))
    train_sampler = torch.utils.data.distributed.DistributedSampler(dataset=train_dataset)

    valid_dataset = datasets.MNIST(DIR, train=False, transform=transforms.ToTensor())
    valid_dataset = torch.utils.data.Subset(valid_dataset, indices=range(N_VALID_EXAMPLES))
    valid_sampler = torch.utils.data.distributed.DistributedSampler(
        dataset=valid_dataset, shuffle=False
    )

    train_loader = torch.utils.data.DataLoader(
        train_dataset,
        sampler=train_sampler,
        batch_size=BATCHSIZE,
        shuffle=False,
    )
    valid_loader = torch.utils.data.DataLoader(
        valid_dataset,
        sampler=valid_sampler,
        batch_size=BATCHSIZE,
        shuffle=False,
    )

    return train_loader, valid_loader, train_sampler, valid_sampler


def objective(single_trial):
    print("starting objective")
    trial = optuna.integration.TorchDistributedTrial(single_trial, rank)

    # Generate the model.
    model = DDP(define_model(trial).to(rank), device_ids=[rank])

    # Generate the optimizers.
    optimizer_name = trial.suggest_categorical("optimizer", ["Adam", "RMSprop", "SGD"])
    lr = trial.suggest_float("lr", 1e-5, 1e-1, log=True)
    optimizer = getattr(optim, optimizer_name)(model.parameters(), lr=lr)

    # Get the MNIST dataset.
    # train_loader, valid_loader, train_sampler, valid_sampler = get_mnist()

    accuracy = 0
    # Training of the model.
    print("starting to train")
    for epoch in range(EPOCHS):
        dist.barrier()
        accuracy = 0
        trial.report(accuracy, epoch)

        # Handle pruning based on the intermediate value.
        if trial.should_prune():
            raise optuna.exceptions.TrialPruned()

    return accuracy


if __name__ == "__main__":
    print("starting")
    parser = argparse.ArgumentParser()
    parser.add_argument("--local_rank", type=int)
    args = parser.parse_args()

    # fix seed
    seed = 7
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)

    # or 2
    world_size = 2

    rank = args.local_rank
    print(rank)

    rank = int(os.environ["LOCAL_RANK"])
    dist.init_process_group(
        backend="nccl",
        init_method="env://"
    )
    device_rank = dist.get_rank()
    torch.cuda.set_device(rank)
    print("ddp set up")
    print("past barrier")

    study = None
    n_trials = 20
    if device_rank == 0:
        study = optuna.create_study(direction="maximize")
        study.optimize(objective, n_trials=n_trials)
    else:
        for i in range(n_trials):
            try:
                print(f"starting trial: {i}")
                objective(None)
            except optuna.TrialPruned:
                pass

    if device_rank == 0:
        assert study is not None
        pruned_trials = [t for t in study.trials if t.state == optuna.trial.TrialState.PRUNED]
        complete_trials = [t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE]

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