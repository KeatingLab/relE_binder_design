#!/bin/bash
  
LLSUB_SIZE=1
LLSUB_RANK=0

echo LLSUB_SIZE $LLSUB_SIZE
echo LLSUB_RANK $LLSUB_RANK

#bridgeDB=/data1/groups/keatinglab/swans/databases/singleChainDB_max8res/job/singlechain_190122_22188.bridge.bin
bridgeDB=/orcd/pool/008/swans/databases/singleChainDB_max8res/job/singlechain_190122_22188.bridge.bin
#structureDB=/data1/groups/keatinglab/swans/databases/singlechain_190122_22188.db
structureDB=/orcd/pool/008/swans/databases/singlechain_190122_22188.db
dCut=1.0
RMSDCut=0.75
seedAMultiPDB=../files/4FXE-ARG81-relax-noHyd-48-66__0.pdb
seedBMultiPDB=../01_generateSeeds/4FXE-ARG81-relax-noHyd_E_B.seeds.pdb
resectLength=1
minSeedLength=5
maxBridgeLength=4
avoidClashesToStructure=../files/4FXE-ARG81-relax-noHyd_E_.pdb
workerID=$LLSUB_RANK
nWorkers=$LLSUB_SIZE

bin=../../fragment_tools_cpp/bin/connectFragmentsMapper

SECONDS=0

echo $bin --bridgeDB $bridgeDB --structureDB $structureDB --seedAMultiPDB $seedAMultiPDB --seedBMultiPDB $seedBMultiPDB --resectLength $resectLength --minSeedLength $minSeedLength --maxBridgeLength $maxBridgeLength --dCut $dCut --RMSDCut $RMSDCut --avoidClashesToStructure $avoidClashesToStructure --workerID $workerID --nWorkers $nWorkers
$bin --bridgeDB $bridgeDB --structureDB $structureDB --seedAMultiPDB $seedAMultiPDB --seedBMultiPDB $seedBMultiPDB --resectLength $resectLength --minSeedLength $minSeedLength --maxBridgeLength $maxBridgeLength --dCut $dCut --RMSDCut $RMSDCut --avoidClashesToStructure $avoidClashesToStructure --workerID $workerID --nWorkers $nWorkers

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED