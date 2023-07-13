#!/bin/bash

source $(conda info --base)/etc/profile.d/conda.sh
conda activate dl_binder_design

relax_path=/data1/groups/keatinglab/swans/binderDesign_relE/analysis/round_4_thesis/rosetta_analysis/1_relax_denovo2seedMPNN
name=denovo_2seed_mpnn
tagslist=/data1/groups/keatinglab/swans/binderDesign_relE/analysis/combined_analysis/filtereddesigns_tags/${name}_list.txt
batch_size=315
    
cat ${relax_path}/*.silent > ${name}.silent

cat $tagslist | silentfasterslice ${name}.silent > ${name}_filtered.silent

mkdir silentsplit
cd silentsplit

silentsplit ../${name}_filtered.silent $batch_size

#rename
files=(*)
total=${#files[@]}
i=0
for f in "${files[@]}"; do
    echo index $i
    echo "- Processing file: $f"
    mv $f ../${name}_$i.silent
    i=$(( i + 1 ))
done

cd ..
rm -r silentsplit
echo "Done!"