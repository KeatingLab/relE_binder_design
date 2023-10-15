#!/bin/bash
. /etc/profile.d/modules.sh

source $(conda info --base)/etc/profile.d/conda.sh
conda activate dl_binder_design

# LLSUB_SIZE=192
# LLSUB_RANK=0

#######################################################################
# Set these variables before running
mpnn_dir=/data1/groups/keatinglab/swans/binderDesign_relE/data/round_4_3seeds/6_designSequenceMPNN_3seeddenovo_5seqperbb/
multientry_pdb=/data1/groups/keatinglab/swans/binderDesign_relE/analysis/round_4_3seeds/multientrypdb_4FXE_3seeds_structscoreperres-.5_seqstructtop5k_modifiedPDBs/4FXE_3seeds_structscoreperres-.5_seqstructtop5k_modifiedPDBs.pdb
silent_name=denovo3seed_mpnn_top5seq
target_pdb=/data1/groups/keatinglab/swans/binderDesign_relE/files/4FXE-ARG81-relax-noHyd_E_.pdb

runfile=/home/gridsan/sswanson/local_code_mirror/peptide_binder_design/python/scripts/createSilentFileMPNN.py
#######################################################################

SECONDS=0

python $runfile \
        --mpnn_dir $mpnn_dir \
        --multientry_pdb $multientry_pdb \
        --silent_name $silent_name \
        --target_pdb $target_pdb \
        --n_batches $LLSUB_SIZE \
        --batch_id $LLSUB_RANK

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED