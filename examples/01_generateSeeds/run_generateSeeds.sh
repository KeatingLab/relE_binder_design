#!/bin/bash 
#SBATCH -J generateSeeds
#SBATCH -n 1 #Request 1 task (core)
#SBATCH -N 1 #Request 1 node
#SBATCH -t 0-3:00
#SBATCH --mem-per-cpu=10G #Request 4G of memory per CPU
#SBATCH -o generateSeeds.%J.out #redirect output to output_JOBID.txt
#SBATCH -e generateSeeds.%J.err #redirect errors to error_JOBID.txt

fasstDB=/orcd/pool/008/swans/databases/220211_PDB_biounits_2019-01-22/nrBioUnits_190122_23633.db
complexPDB=../files/4FXE-ARG81-relax-noHyd_E_B.pdb
binderChains="B"
numMatches=10 #small for the example. Originally used 100000
seqConst=""
pdbOut="--saveToPDB"
seedFlankRes=3
wholeSurface="--wholeSurface"

SECONDS=0

bin=../../fragment_tools_cpp/bin/generateSeeds

echo $bin --complexPDB $complexPDB --binderChains $binderChains --fasstDB $fasstDB --numMatches $numMatches --seedFlankRes $seedFlankRes $seqConst $pdbOut $wholeSurface
$bin --complexPDB $complexPDB --binderChains $binderChains --fasstDB $fasstDB --numMatches $numMatches --seedFlankRes $seedFlankRes $seqConst $pdbOut $wholeSurface

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED
