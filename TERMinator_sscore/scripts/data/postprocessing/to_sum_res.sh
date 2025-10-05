#!/bin/bash
#SBATCH -N 1
#SBATCH --partition=xeon-p8
#SBATCH --time=30:00:00
#SBATCH --mem=0
#SBATCH -o OUTPUTDIR/res-output.out
#SBATCH -e OUTPUTDIR/res-error.out

# activate conda
CONDA_ROOT=/state/partition1/llgrid/pkg/anaconda/anaconda3-2022b/
source ${CONDA_ROOT}/etc/profile.d/conda.sh
module load anaconda/2022b
conda activate terminator-nightly
echo OUTPUTDIR

python summarize_results.py \
    --output_dir=OUTPUTDIR \
    --dtermen_data=DTERMENDATA \
