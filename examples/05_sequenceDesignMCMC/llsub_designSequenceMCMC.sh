#!/bin/bash

source $(conda info --base)/etc/profile.d/conda.sh
conda activate terminator_sscore

export PYTHONPATH=$PYTHONPATH:~/local_code_mirror/peptide_binder_design
export PYTHONPATH=$PYTHONPATH:~/local_code_mirror/TERMinator_sscore
repo=~/local_code_mirror/peptide_binder_design

#LLSUB_RANK=0
#LLSUB_SIZE=10
echo $LLSUB_RANK
echo $LLSUB_SIZE

SECONDS=0

python -u $repo/python/design_sequence/designSequenceMCMC.py \
      --complex_dataset $PWD/bridgeSeeds_5D94_100kseeds_top1000totalscore_wtarget_modifiedPDBs_list.txt \
      --binderChainID "A" \
      --targetChainID "Z" \
      --model_dir /data1/groups/keatinglab/swans/TERMinator/finetuneCOORDinator_sscore_experiments/selfedge_MLP_64 \
      --early_stopping 3 \
      --n 1 \
      --n_batches $LLSUB_SIZE \
      --batch_index $LLSUB_RANK

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED