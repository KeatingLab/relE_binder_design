#!/bin/bash
# SBATCH -N 1
# SBATCH --job-name clusterBindersByContacts
# SBATCH --tasks-per-node=1
# SBATCH --mem-per-cpu=4G #Request 4G of memory per CPU
# SBATCH -e clusterBindersByContacts.%j.err
# SBATCH -o clusterBindersByContacts.%j.out

contactData=/home/gridsan/sswanson/local_code_mirror/peptide_binder_design/testfiles/contact2DHistogramData_withClashData_withUNK.json
targetPDB=/data1/groups/keatinglab/swans/binderDesign_relE/files/4FXE-ARG81-relax-noHyd_E_.pdb
structureList=4FXE_3seeds_structscoreperres-.5_seqstructtop5k_modifiedPDBs_list.txt
#structureList=../4FXE_denovoextensions_10k_top5k-structscore-.25-seqstructscoreperres_MCMC_wtarget_top10_list.txt
binderChainID="0"
cutoff=0.25
binderChainResidueDistance=12
batchID=$LLSUB_RANK
nBatches=$LLSUB_SIZE

bin=/home/gridsan/sswanson/local_code_mirror/peptide_binder_design/cplusplus/bin/clusterBindersByContacts

SECONDS=0

echo $bin --contactData $contactData --targetPDB $targetPDB --structureList $structureList --binderChainID $binderChainID --cutoff $cutoff --binderChainResidueDistance $binderChainResidueDistance --batchID $batchID --nBatches $nBatches
$bin --contactData $contactData --targetPDB $targetPDB --structureList $structureList --binderChainID $binderChainID --cutoff $cutoff --binderChainResidueDistance $binderChainResidueDistance --batchID $batchID --nBatches $nBatches

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED