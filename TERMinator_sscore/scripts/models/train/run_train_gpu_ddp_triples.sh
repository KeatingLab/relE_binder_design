#!/bin/bash

# set up conda environment
CONDA_ROOT=/state/partition1/llgrid/pkg/anaconda/anaconda3-2022b/
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate terminator_sscore
echo $CONDA_PREFIX

# set variables
NODERANK=$1
MASTER_ADDR=$2
DATASET=$3
MODEL_HPARAMS=$4
RUN_HPARAMS=$5
RUNDIR=$6
TERMINATOR_REPO=$7
TRAIN=$DATASET/$8
VALIDATION=$DATASET/$9
TEST=$DATASET/${10}
MASTER_PORT=1234

echo $NODERANK
echo $TRAIN
echo $VALIDATION
echo $TEST

# set modules and core dump params
ulimit -s unlimited
ulimit -n 10000
module load cuda/11.1
module load nccl/2.8.3-cuda11.1
export NCCL_DEBUG=INFO

train_script=$TERMINATOR_REPO/scripts/models/train/train.py

python -m torch.distributed.run --nnodes $NUM_NODES --nproc_per_node $NUM_GPUS_PER_NODE --node_rank $NODERANK --master_addr $MASTER_ADDR --master_port "$MASTER_PORT" --rdzv_id 2 --rdzv_backend c10d --rdzv_endpoint $MASTER_ADDR:"$MASTER_PORT" $train_script --dataset=$DATASET --model_hparams=$MODEL_HPARAMS --backend=nccl --run_hparams=$RUN_HPARAMS --run_dir=$RUNDIR --n_nodes=$NUM_NODES --n_trials=$NUM_TRIALS --train $TRAIN --validation $VALIDATION --test $TEST


