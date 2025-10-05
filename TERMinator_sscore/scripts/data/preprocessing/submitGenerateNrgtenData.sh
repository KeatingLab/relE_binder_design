#!/bin/bash
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --partition=xeon-p8
#SBATCH --time=30:00:00
#SBATCH -o ./dynamic-output-valid.out
#SBATCH -e ./dynamic-error-valid.out

CONDA_ROOT=/state/partition1/llgrid/pkg/anaconda/anaconda3-2022b/
source ${CONDA_ROOT}/etc/profile.d/conda.sh
module load anaconda/2022b
conda activate terminator-nightly

python ~/TERMinator/scripts/data/preprocessing/generateNrgtenData.py --pdb_root /data1/groups/keatinglab/fosterb/data/multichain_clean --pdb_list /data1/groups/keatinglab/fosterb/data/multichain_nrgten_small_motion_78/validation.in --out_dir ~/keatinglab_shared/fosterb/data/multichain_nrgten_large_motion_78 --nrgten_type signatures_ensembles --num_processes 64 --success_file validation.in --update False