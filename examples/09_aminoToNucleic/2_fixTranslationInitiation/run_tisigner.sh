#!/bin/bash
. /etc/profile.d/modules.sh

#NOTE: must add the following line to ~/.bashrc for this script to work
#PATH=$PATH:/data1/groups/keatinglab/swans/repos/ViennaRNA/bin
#This allows the script to call the external application for calculating probabilities of RNA structures

################################################################
# parameters to set
utr="tccctatcagtgatagagattgacatccctatcagtgatagagatactgagcacctaggtttctccataccttatcaaacgggactcaaatctt"
host="ecoli"

# dna sequence (from dnaworks, should include ATG)
sequence=$1

# name of output file
output=$2

################################################################

source $(conda info --base)/etc/profile.d/conda.sh
conda activate tisigner

# parameters that shouldn't be changed
codons=9
seed=42

prog=/data1/groups/keatinglab/swans/repos/TIsigner/TIsigner_cmd/tisigner.py

SECONDS=0

python $prog -s $sequence -o $output -c $codons -u $utr -h $host -d $seed

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED