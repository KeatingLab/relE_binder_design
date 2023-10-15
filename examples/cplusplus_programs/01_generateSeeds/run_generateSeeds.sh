#!/bin/bash 
#SBATCH -J generateSeeds
#SBATCH -n 1 #Request 1 task (core)
#SBATCH -N 1 #Request 1 node
#SBATCH -t 0-3:00
#SBATCH -C centos7 #Request only Centos7 nodes
#SBATCH -p sched_mit_hill #Run on sched_engaging_default partition
#SBATCH --mem-per-cpu=16G #Request 4G of memory per CPU
#SBATCH -o generateSeeds%j.out #redirect output to output_JOBID.txt
#SBATCH -e generateSeeds%j.err #redirect errors to error_JOBID.txt

fasstDB=/nobackup1/swans/resFrameScoring/buildDB/220211_PDB_biounits_2019-01-22/nrBioUnits_190122_23633.db
complexPDB=../../testfiles/5UUL.pdb
binderChains="B"
numMatches=10000
seqConst=""

SECONDS=0

echo /home/swans/MST_repos/interfaceGenerator/bin/generateSeeds --complexPDB $complexPDB --binderChains $binderChains --fasstDB $fasstDB --numMatches $numMatches $seqConst
/home/swans/MST_repos/interfaceGenerator/bin/generateSeeds --complexPDB $complexPDB --binderChains $binderChains --fasstDB $fasstDB --numMatches $numMatches $seqConst

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED
