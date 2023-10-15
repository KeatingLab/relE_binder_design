#!/bin/bash

# LLsub ./llsub_designSequenceMCMC.sh [1,48,1] -q xeon-p8

source $(conda info --base)/etc/profile.d/conda.sh
conda activate terminator_sscore

export PYTHONPATH=$PYTHONPATH:~/local_code_mirror/peptide_binder_design
export PYTHONPATH=$PYTHONPATH:~/local_code_mirror/TERMinator_sscore

repo=~/local_code_mirror/peptide_binder_design

#LLSUB_RANK=0
#LLSUB_SIZE=1
echo $LLSUB_RANK
echo $LLSUB_SIZE

SECONDS=0

python -u $repo/python/design_sequence/designSequenceMCMC.py \
      --complex_dataset 4FXE_3seeds_structscoreperres-.5_seqstructtop5k_modifiedPDBs_list.txt \
      --binderChainID "0" \
      --targetChainID "E" \
      --model_dir /data1/groups/keatinglab/swans/TERMinator/230310_finetuneCOORDinator_sscore_mpnodeupdate \
      --early_stopping -1 \
      --n 1 \
      --n_cyc 100 \
      --n_it 10000 \
      --Tf 0.5 \
      --n_batches $LLSUB_SIZE \
      --batch_index $LLSUB_RANK

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED