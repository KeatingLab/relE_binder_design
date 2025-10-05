#!/bin/bash
#SBATCH -N 1
#SBATCH --mincpu=32
#SBATCH --gres=gpu:volta:1
#SBATCH --time=HOURS:00:00
#SBATCH --mem=60G
#SBATCH -o RUNDIR/train-output_runRUNNO.out
#SBATCH -e RUNDIR/train-error_runRUNNO.out

CONDA_ROOT=/state/partition1/llgrid/pkg/anaconda/anaconda3-2022b/
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate terminator_sscore
ulimit -s unlimited
ulimit -n 10000

echo $CONDA_PREFIX

python train.py \
  --dataset=DATASET \
  --model_hparams=MODEL_HPARAMS \
  --run_hparams=RUN_HPARAMS \
  --run_dir=RUNDIR \
  --out_dir=OUTPUTDIR \
  --lazy
