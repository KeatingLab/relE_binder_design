#!/bin/bash
#SBATCH -N 1
#SBATCH --partition=xeon-p8
#SBATCH --time=30:00:00
#SBATCH -o ./convert-etab-out_e.out
#SBATCH -e ./convert-etab-error_e.out

CONDA_ROOT=/state/partition1/llgrid/pkg/anaconda/anaconda3-2022b/
ROOT_FOLDER=~/TERMinator_experiments/flex_ensembles_embedded/rocklin_results
source ${CONDA_ROOT}/etc/profile.d/conda.sh
module load anaconda/2022b
conda activate terminator-nightly

python ~/TERMinator/scripts/data/postprocessing/condenseEtabs.py --in_folder ${ROOT_FOLDER}/etabs --out_folder ${ROOT_FOLDER}/etabs_condensed/ -n 64