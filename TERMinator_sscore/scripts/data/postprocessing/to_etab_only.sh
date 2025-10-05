#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --partition=xeon-p8
#SBATCH --time=20:00:00
#SBATCH --mem=0
#SBATCH -o OUTPUTDIR/etab-output.out
#SBATCH -e OUTPUTDIR/etab-error.out
#SBATCH --exclusive

# activate conda
CONDA_ROOT=/state/partition1/llgrid/pkg/anaconda/anaconda3-2022b/
source ${CONDA_ROOT}/etc/profile.d/conda.sh
module load anaconda/2022b
conda activate terminator-nightly

python to_etab.py \
    --output_dir=OUTPUTDIR \
    --dtermen_data=DTERMENDATA \
    --num_cores=64 -u
