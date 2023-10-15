#!/bin/bash
#SBATCH -J clusterSeeds
#SBATCH -n 1 #Request 1 task (core)
#SBATCH -N 1 #Request 1 node
#SBATCH -t 0-3:00
#SBATCH --mem-per-cpu=10G #Request 4G of memory per CPU
#SBATCH -o clusterSeeds.%J.out #redirect output to output_JOBID.txt
#SBATCH -e clusterSeeds.%J.err #redirect errors to error_JOBID.txt

bin=/home/gridsan/sswanson/local_code_mirror/peptide_binder_design/cplusplus/bin/greedyClusterSeeds

seedBin=/data1/groups/keatinglab/swans/binderDesign_lc3b/round_1/data/1_generateSeeds_5D94_wholeSurface_100k_fixed/5D94_A_B.seeds.bin
# targetPDB=/data1/groups/keatinglab/swans/binderDesign_lc3b/round_1/files/5D94_A__.pdb
seedList=/data1/groups/keatinglab/swans/binderDesign_lc3b/round_2/analysis/230304_top10kseeds_totalscore_list.txt
clusterRMSDCutoff=1.0
seedLength=7
coverage=.75

SECONDS=0

echo $bin --binderBin $seedBin --seedList $seedList --clusterRMSDCutoff $clusterRMSDCutoff --seedLength $seedLength --coverage $coverage
$bin --binderBin $seedBin --seedList $seedList --clusterRMSDCutoff $clusterRMSDCutoff --seedLength $seedLength --coverage $coverage

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED