#!/bin/bash
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --partition=xeon-p8
#SBATCH --time=30:00:00
#SBATCH --mem=0
#SBATCH -o OUTPUTDIR/etab-output.out
#SBATCH -e OUTPUTDIR/etab-error.out

# activate conda
CONDA_ROOT=/state/partition1/llgrid/pkg/anaconda/anaconda3-2022b/
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate terminator

python batch_arr_dTERMen.py \
    --output_dir=OUTPUTDIR \
    --pdb_root=PDBROOT \
    --dtermen_data=DTERMENDATA \
    --batch_size=48