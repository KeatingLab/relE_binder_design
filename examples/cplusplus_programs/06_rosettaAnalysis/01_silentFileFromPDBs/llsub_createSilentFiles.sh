#!/bin/bash
. /etc/profile.d/modules.sh

source $(conda info --base)/etc/profile.d/conda.sh
conda activate dl_binder_design
#export PYTHONPATH=$PYTHONPATH:~/local_code_mirror/peptide_binder_design
#export PYTHONPATH=$PYTHONPATH:~/local_code_mirror/TERMinator_sscore

#LLSUB_SIZE=192
#LLSUB_RANK=0

#######################################################################
# Set these variables before running
selected_designs=/data1/groups/keatinglab/swans/binderDesign_relE/analysis/round_4_thesis/230510_denovo_top5.csv
multientry_pdb=/data1/groups/keatinglab/swans/binderDesign_relE/data/round_4/5_scoreSeeds_bridgeSeeds_top100dGSepSeeds_top10kseeds_proteomebg/nterm_cterm_10kseeds_fused.pdb
silent_name=denovo_top5seq
target_pdb=/data1/groups/keatinglab/swans/binderDesign_relE/files/4FXE-ARG81-relax-noHyd_E_.pdb
binder_chain_id="0"
name_col=bridge_name
seq_col=mcmc_seq_peptide

runfile=/home/gridsan/sswanson/local_code_mirror/peptide_binder_design/python/pdb_processing/createSilentFile.py
#######################################################################

python $runfile \
	--selected_designs $selected_designs \
 	--multientry_pdb $multientry_pdb \
  	--silent_name $silent_name \
  	--target_pdb $target_pdb \
	--n_batches $LLSUB_SIZE \
	--batch_id $LLSUB_RANK \
	--binder_chain_id $binder_chain_id \
	--name_col $name_col \
	--seq_col $seq_col
