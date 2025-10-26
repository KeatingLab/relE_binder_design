#!/bin/bash
. /etc/profile.d/modules.sh

source $(conda info --base)/etc/profile.d/conda.sh
conda activate dl_binder_design

#######################################################################
# Set these variables before running
mpnn_dir=6_designSequenceMPNN_3seeddenovo_5seqperbb # need to run MPNN first
multientry_pdb=4FXE_3seeds_structscoreperres-.5_seqstructtop5k_modifiedPDBs.pdb
silent_name=denovo3seed_mpnn_top5seq
target_pdb=../files/4FXE-ARG81-relax-noHyd_E_.pdb

runfile=../../scripts/python/createSilentFileMPNN.py
#######################################################################

SECONDS=0

python $runfile \
        --mpnn_dir $mpnn_dir \
        --multientry_pdb $multientry_pdb \
        --silent_name $silent_name \
        --target_pdb $target_pdb \
        --n_batches 1 \
        --batch_id 0

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED