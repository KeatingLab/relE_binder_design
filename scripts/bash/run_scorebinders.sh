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
conda activate terminator_sscore
export PYTHONPATH=$PYTHONPATH:~/local_code_mirror/peptide_binder_design
export PYTHONPATH=$PYTHONPATH:~/local_code_mirror/TERMinator_sscore

binder_dataset=
TERMrepo=/home/gridsan/sswanson/local_code_mirror/peptide_binder_design
model=/data1/groups/keatinglab/swans/TERMinator/230218_finetuneCOORDinator_sscore_selfedge
seqMode="consensus_aa"
scoreMode="interface_only"

SECONDS=0

python $TERMrepo/scripts/score_binders/scoreBinders.py \
      --binder_dataset $binder_dataset \
      --model_dir $model \
      --seq_mode $seqMode \
      --score_mode $scoreMode 

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED
