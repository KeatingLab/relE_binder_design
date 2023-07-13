#!/bin/bash
. /etc/profile.d/modules.sh

source $(conda info --base)/etc/profile.d/conda.sh
conda activate dl_binder_design

#LLSUB_SIZE=192
#LLSUB_RANK=0

echo $LLSUB_SIZE
echo $LLSUB_RANK

#######################################################################
# Set these variables before running
silent_file=/data1/groups/keatinglab/swans/binderDesign_relE/analysis/round_4_thesis/rosetta_analysis/0_silentFromPDBs_denovo2seedMPNN/denovo2seed_mpnn_top5seq_${LLSUB_RANK}_binder_designs.silent
target_sequence="AYFLDFDERALKEWRKLGSTVREQLKKKLVEVLESPRIEANKLRGMPDCYKIKLRSSGYRLVYQVIDEKVVVFVISVGKRER"

runfile=/home/gridsan/sswanson/local_code_mirror/peptide_binder_design/python/scripts/silentFileExtractSequence.py
#######################################################################

python $runfile \
	--silent_file $silent_file \
    --target_sequence "$target_sequence" \
    --batch_id $LLSUB_RANK