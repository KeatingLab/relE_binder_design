#!/bin/bash
#SBATCH -N 1
#SBATCH --mincpu=32
#SBATCH --gres=gpu:volta:1
#SBATCH --time=10:00:00
#SBATCH --mem=60G
#SBATCH -o ~/train-output_run_test.out
#SBATCH -e ~/train-error_run_test.out

CONDA_ROOT=/state/partition1/llgrid/pkg/anaconda/anaconda3-2022b/
source ${CONDA_ROOT}/etc/profile.d/conda.sh

time python -m torch.distributed.launch --nproc_per_node=2 pytorch_distributed_simply.py
