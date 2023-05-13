#!/bin/bash
#SBATCH -J seedBridge
#SBATCH -n 1 #Request 1 task (core)
#SBATCH -N 1 #Request 1 node
#SBATCH -t 1-0:00:00
#SBATCH --mem-per-cpu=16G #Request 4G of memory per CPU
#SBATCH -o seedBridge.%J.out #redirect output to output_JOBID.txt
#SBATCH -e seedBridge.%J.err #redirect errors to error_JOBID.txt

bridgeDB=/home/gridsan/sswanson/keatinglab_shared/swans/databases/singleChainDB_max8res/job/singlechain_190122_22188.bridge.bin
structureDB=/data1/groups/keatinglab/swans/databases/singlechain_190122_22188.db
dCut=1.0
RMSDCut=0.75
seedBin=/data1/groups/keatinglab/swans/binderDesign_lc3b/round_1/data/1_generateSeeds_5D94_wholeSurface_100k_fixed/5D94_A_B.seeds.bin
seedList=/data1/groups/keatinglab/swans/binderDesign_lc3b/round_2/analysis/230304_top10kseeds_totalscore_list.txt
seedA=/data1/groups/keatinglab/swans/binderDesign_lc3b/round_2/files/5D94__B.pdb
resectLength=4
minSeedLength=4
avoidClashesToStructure=/data1/groups/keatinglab/swans/binderDesign_lc3b/round_2/files/5D94_A_.pdb

bin=/home/gridsan/sswanson/local_code_mirror/peptide_binder_design/cplusplus/bin/buildBridgeSeedsGraph

SECONDS=0

echo $bin $bridgeDB --structureDB $structureDB --seedBin $seedBin --seedList $seedList --seedA $seedA --resectLength $resectLength --minSeedLength $minSeedLength --dCut $dCut --RMSDCut $RMSDCut --avoidClashesToStructure $avoidClashesToStructure
$bin --bridgeDB $bridgeDB --structureDB $structureDB --seedBin $seedBin --seedList $seedList --seedA $seedA --resectLength $resectLength --minSeedLength $minSeedLength --dCut $dCut --RMSDCut $RMSDCut --avoidClashesToStructure $avoidClashesToStructure

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED