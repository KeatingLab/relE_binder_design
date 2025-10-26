#!/bin/bash

export PYTHONPATH=$PYTHONPATH:/orcd/home/002/swans/relE_paper/relE_binder_design/
export PYTHONPATH=$PYTHONPATH:/orcd/home/002/swans/relE_paper/relE_binder_design/TERMinator_sscore

repo=../..

complex_dataset=relBE.txt
model=/orcd/pool/008/swans/terminator_stuff/priority/230310_finetuneCOORDinator_sscore_mpnodeupdate

SECONDS=0

poetry run python -u $repo/scripts/python/designSequenceMCMC.py \
      --complex_dataset $complex_dataset \
      --binderChainID "B" \
      --targetChainID "E" \
      --model_dir $model \
      --early_stopping -1 \
      --n 1 \
      --n_cyc 5 \
      --n_it 10000 \
      --Tf 0.5 \
      --n_batches 1 \
      --batch_index 0

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED
