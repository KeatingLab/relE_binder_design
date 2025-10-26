#!/bin/bash
#SBATCH -J rosettaRelax
#SBATCH -n 1 #Request 1 task (core)
#SBATCH -N 1 #Request 1 node
#SBATCH -t 0-3:00
#SBATCH --mem-per-cpu=4G #Request 4G of memory per CPU
#SBATCH -o rosettaRelax.%J.out #redirect output to output_JOBID.txt
#SBATCH -e rosettaRelax.%J.err #redirect errors to error_JOBID.txt

silentfile=$1
silentfile_name=${silentfile##*/}
relaxScript=../../scripts/rosetta/relaxscript_modified_InterfaceRelax2019_constrainttostart.txt
args=""

rosettaDir=/data1/groups/keatinglab/rosetta/rosetta_src_2021.16.61629_bundle/main
relaxBin=$rosettaDir/source/bin/relax.linuxgccrelease

SECONDS=0
echo $relaxBin -relax:script $relaxScript -in:file:silent $silentfile -out:file:silent $silentfile_name -out:file:scorefile score.sc $args
$relaxBin -relax:script $relaxScript -in:file:silent $silentfile -out:file:silent $silentfile_name -out:file:scorefile score.sc $args
ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED
