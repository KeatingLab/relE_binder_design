#!/bin/bash
#SBATCH -N 1
#SBATCH --time=01:00:00
#SBATCH -J submit_design
#SBATCH -o submit_design.%j.out
#SBATCH -e submit_design.%j.err

path=/home/gridan/fbirnbaum/TERMinator_experiments/noise_batch_0/run_noise_level_0.0/results/etabs
sbatch=/home/gridan/fbirnbaum/TERMinator/scripts/data/postprocessing/run_designsequence.sh

for etab in $(ls $path/*.etab); do
	echo $etab
	file=${etab##*/}
	name=${file%.etab}
	mkdir $name
	cd $name
	sbatch $sbatch $etab
	cd ..
done
