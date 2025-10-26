#!/bin/bash
# SBATCH -N 1
# SBATCH --job-name clusterBindersByContacts
# SBATCH --tasks-per-node=1
# SBATCH --mem-per-cpu=4G #Request 4G of memory per CPU
# SBATCH -e clusterBindersByContacts.%j.err
# SBATCH -o clusterBindersByContacts.%j.out

contactData=/home/gridsan/sswanson/local_code_mirror/peptide_binder_design/testfiles/contact2DHistogramData_withClashData_withUNK.json
targetPDB=/data1/groups/keatinglab/swans/binderDesign_relE/files/4FXE-ARG81-relax-noHyd_E_.pdb
multiPDB=/data1/groups/keatinglab/swans/binderDesign_relE/analysis/round_4_2seedextensions/relBhelix_1seedextensions.pdb
binderChainID="0"
cutoff=0.25
binderChainResidueDistance=12
batchID=0
nBatches=1

bin=../../fragment_tools_cpp/bin/clusterBindersByContacts

SECONDS=0

echo $bin --contactData $contactData --targetPDB $targetPDB --multiPDB $multiPDB --binderChainID $binderChainID --cutoff $cutoff --binderChainResidueDistance $binderChainResidueDistance --batchID $batchID --nBatches $nBatches
$bin --contactData $contactData --targetPDB $targetPDB --multiPDB $multiPDB --binderChainID $binderChainID --cutoff $cutoff --binderChainResidueDistance $binderChainResidueDistance --batchID $batchID --nBatches $nBatches

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED
