#!/bin/bash
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --partition=xeon-p8
#SBATCH --time=30:00:00
#SBATCH -o ./nrgten-validate-out-slb.out
#SBATCH -e ./nrgten-validate-error-slb.out

CONDA_ROOT=/state/partition1/llgrid/pkg/anaconda/anaconda3-2022b/
source ${CONDA_ROOT}/etc/profile.d/conda.sh
module load anaconda/2022b
conda activate terminator-nightly

python ~/TERMinator/scripts/data/preprocessing/validateNrgtenData.py --full_root ~/keatinglab_shared/fosterb/data/multichain_nrgten_large_motion_789 --full_list ~/keatinglab_shared/fosterb/data/multichain_nrgten_large_motion_789/validation.in --backbone_root ~/keatinglab_shared/fosterb/data/multichain_nrgten_small_motion_78 --backbone_list ~/keatinglab_shared/fosterb/data/multichain_nrgten_small_motion_78/validation.in --out_dir ~/nrgten_examples/comp_results_small_large --pdb_dir ~/keatinglab_shared/fosterb/data/multichain_clean --num_processes 64 --analysis_level atom --verbose True --update True --analysis_type ensembles_b_factors