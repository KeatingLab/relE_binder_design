#!/bin/bash
#SBATCH -N 1
#SBATCH --mincpu=32
#SBATCH --gres=gpu:volta:1
#SBATCH --time=HOURS:00:00
#SBATCH --mem=60G
#SBATCH -o RUNDIR/train-output_runRUNNO_NODERANK.out
#SBATCH -e RUNDIR/train-error_runRUNNO_NODERANK.out

CONDA_ROOT=/state/partition1/llgrid/pkg/anaconda/anaconda3-2022b/
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate terminator_sscore

MASTER_PORT=1234
ulimit -s unlimited
ulimit -n 10000
echo NODERANK
module load cuda/11.1
module load nccl/2.8.3-cuda11.1
conda activate terminator-nightly
export NCCL_DEBUG=INFO

echo $CONDA_PREFIX

python -m torch.distributed.run --nnodes NUM_NODES --nproc_per_node NUM_GPUS_PER_NODE --node_rank NODERANK --master_addr MASTER_ADDR --master_port "$MASTER_PORT" --rdzv_id 2 --rdzv_backend c10d --rdzv_endpoint MASTER_ADDR:"$MASTER_PORT" train.py --dataset=DATASET --model_hparams=MODEL_HPARAMS --backend=nccl --run_hparams=RUN_HPARAMS --run_dir=RUNDIR --n_nodes=NUM_NODES --n_trials=NUM_TRIALS --train=DATASET/TRAIN --validation=DATASET/VALIDATION --test=DATASET/TEST


