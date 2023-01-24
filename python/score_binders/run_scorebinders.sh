#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name scoreBinders
#SBATCH --tasks-per-node=1
#SBATCH --mem-per-cpu=8G #Request 4G of memory per CPU
#SBATCH -e scoreBinders.%j.err
#SBATCH -o scoreBinders.%j.out
#SBATCH -t 0-04:00:00
#SBATCH --gres gpu:volta:1

source $(conda info --base)/etc/profile.d/conda.sh
conda activate terminator
export PYTHONPATH=~/local_code_mirror/TERMinator

binder_dataset=
TERMrepo=/home/gridsan/sswanson/local_code_mirror/TERMinator
model=$TERMrepo/models/multichain_coordinator
scoreMode="consensus_aa"

SECONDS=0

python $TERMrepo/scripts/score_binders/scoreBinders.py \
      --binder_dataset $binder_dataset \
      --model_dir $model \
      --score_mode $scoreMode 

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED
