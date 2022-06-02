#!/bin/bash 
#SBATCH -J scoreBinder
#SBATCH -n 1 #Request 1 task (core)
#SBATCH -N 1 #Request 1 node
#SBATCH -t 0-4:00
#SBATCH -C centos7 #Request only Centos7 nodes
#SBATCH -p sched_mit_hill #Run on sched_engaging_default partition
#SBATCH --mem-per-cpu=10G #Request 4G of memory per CPU
#SBATCH -o scoreBinder.%j.out #redirect output to output_JOBID.txt
#SBATCH -e scoreBinder.%j.err #redirect errors to error_JOBID.txt

# Note: there's a couple different options for providing structures to score
# --binderBin : a seed binary file
# --binderList : a file where each line is the path a binder structure
# --complexPDB : a PDB file containing a target protein and binder structure
# see scoreBinder options for more info

frameDB=/nobackup1/swans/resFrameScoring/buildFrameDB/220314_1.2VDWUB_PDBbiounits_PixelDB30removed_2019-01-22/nrBioUnits_190122_22531_PixelDB30removed.frames.db
contactData=../../testfiles/contact2DHistogramData_withClashData.json
complexPDB=../../testfiles/5UUL.pdb
pepChain="B"
binderBin=../01_generateSeeds/5UUL.seeds.bin

#These parameters seem to work best for scoring based on preliminary tests
distanceCutoff=1.0
orientationCutoff=10

SECONDS=0

echo /home/swans/MST_repos/interfaceGenerator/bin/scoreBinder --frameDB $frameDB --contactData $contactData --targetPDB $targetPDB --binderChains $pepChain  --binderBin $binderBin --distanceCutoff $distanceCutoff --orientationCutoff $orientationCutoff
/home/swans/MST_repos/interfaceGenerator/bin/scoreBinder --frameDB $frameDB --contactData $contactData --targetPDB $targetPDB --binderChains $pepChain --binderBin $binderBin --distanceCutoff $distanceCutoff --orientationCutoff $orientationCutoff

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED
