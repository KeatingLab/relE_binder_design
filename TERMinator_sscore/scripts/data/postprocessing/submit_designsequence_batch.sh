#!/bin/bash
#SBATCH -N 1
#SBATCH --time=01:00:00
#SBATCH -J submit_design
#SBATCH -o submit_design.%j.out
#SBATCH -e submit_design.%j.err

path=/home/gridsan/fbirnbaum/TERMinator_experiments/noise_torsion_batch_0.1_0.1/multichain_results_trim/multichain_results_trim_batch_arr_0.list
etabs=/home/gridsan/fbirnbaum/TERMinator_experiments/noise_torsion_batch_0.1_0.1/multichain_results_trim/etabs
cd $etabs

for pdb in $(cat $path); do
	echo $pdb
	name=run_${pdb}.sh
	echo $name
	sbatch $name
done
