#!/bin/bash
  
#LLSUB_SIZE=1
#LLSUB_RANK=0

echo LLSUB_SIZE $LLSUB_SIZE
echo LLSUB_RANK $LLSUB_RANK

bridgeDB=/data1/groups/keatinglab/swans/databases/singleChainDB_max8res/job/singlechain_190122_22188.bridge.bin
structureDB=/data1/groups/keatinglab/swans/databases/singlechain_190122_22188.db
dCut=1.0
RMSDCut=0.75
seedBin=/data1/groups/keatinglab/swans/binderDesign_relE/data/round_4/1_generateSeeds_100k_noSeqConst_new/4FXE-ARG81-relax-noHyd_E_B.seeds.bin
seedBNames=/data1/groups/keatinglab/swans/binderDesign_relE/analysis/round4/230317_top20kseeds_clusteredto10k_list.txt
seedAPDBList=/data1/groups/keatinglab/swans/binderDesign_relE/files/relBhelix_list.txt
resectLength=1
minSeedLength=5
maxBridgeLength=4
avoidClashesToStructure=/data1/groups/keatinglab/swans/binderDesign_relE/files/4FXE-ARG81-relax-noHyd_E_.pdb
workerID=$LLSUB_RANK
nWorkers=$LLSUB_SIZE

bin=../../fragment_tools_cpp/bin/connectFragmentsMapper

SECONDS=0

echo $bin --bridgeDB $bridgeDB --structureDB $structureDB --seedBin $seedBin --seedBNames $seedBNames --seedAPDBList $seedAPDBList --resectLength $resectLength --minSeedLength $minSeedLength --maxBridgeLength $maxBridgeLength --dCut $dCut --RMSDCut $RMSDCut --avoidClashesToStructure $avoidClashesToStructure --workerID $workerID --nWorkers $nWorkers
$bin --bridgeDB $bridgeDB --structureDB $structureDB --seedBin $seedBin --seedBNames $seedBNames --seedAPDBList $seedAPDBList --resectLength $resectLength --minSeedLength $minSeedLength --maxBridgeLength $maxBridgeLength --dCut $dCut --RMSDCut $RMSDCut --avoidClashesToStructure $avoidClashesToStructure --workerID $workerID --nWorkers $nWorkers

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED