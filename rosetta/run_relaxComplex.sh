#!/bin/bash
#SBATCH -J rosettaRelax
#SBATCH -n 1 #Request 1 task (core)
#SBATCH -N 1 #Request 1 node
#SBATCH -t 0-3:00
#SBATCH --mem-per-cpu=4G #Request 4G of memory per CPU
#SBATCH -o rosettaRelax.%J.out #redirect output to output_JOBID.txt
#SBATCH -e rosettaRelax.%J.err #redirect errors to error_JOBID.txt

inputStructureList=PATH/TO/STRUCTURE
relaxScript=PATH/TO/SCRIPT
args="-relax:constrain_relax_to_native_coords"

rosettaDir=/data1/groups/keatinglab/rosetta/rosetta_src_2021.16.61629_bundle/main
#relaxScript=$rosettaDir/database/sampling/relax_scripts/InterfaceRelax2019.txt
relaxBin=$rosettaDir/source/bin/relax.linuxgccrelease

SECONDS=0
echo $relaxBin -relax:script $relaxScript -in:file:l $inputStructureList $args
$relaxBin -relax:script $relaxScript -in:file:l $inputStructureList $args
ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED