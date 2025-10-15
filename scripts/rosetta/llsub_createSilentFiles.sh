#!/bin/bash
. /etc/profile.d/modules.sh

#######################################################################
# Set these variables before running
select_designs=??? #CSV file containing the names of structures in the multientry PDB (name_col) and peptide sequence (seq_col)
multientry_pdb=/data1/groups/keatinglab/swans/binderDesign_relE/data/round_4/5_scoreSeeds_bridgeSeeds_top100dGSepSeeds_top10kseeds_proteomebg/nterm_cterm_10kseeds_fused.pdb
silent_name=denovo_top5seq
target_pdb=/data1/groups/keatinglab/swans/binderDesign_relE/files/4FXE-ARG81-relax-noHyd_E_.pdbs
binder_chain_id="0"
name_col=bridge_name
seq_col=mcmc_seq_peptide
runfile=$PWD/createSilentFiles.py # the script that will be run for each element in the file
#######################################################################

python $runfile \
  --selected_designs $select_designs \
  --multientry_pdb $multientry_pdb \
  --silent_name $silent_name \
  --target_pdb $target_pdb \
  --n_batches $LLSUB_SIZE \
  --batch_id $LLSUB_RANK \
  --binder_chain_id $binder_chain_id \
  --name_col $name_col \
  --seq_col  $seq_col
