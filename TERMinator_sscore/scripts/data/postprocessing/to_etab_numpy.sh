#!/bin/bash
#SBATCH -N 1
#SBATCH --partition=xeon-p8
#SBATCH --time=30:00:00
#SBATCH --mem=0
#SBATCH -o OUTPUTDIR/etab-nmpy-output.out
#SBATCH -e OUTPUTDIR/etab-nmpy-error.out

# activate conda
CONDA_ROOT=/state/partition1/llgrid/pkg/anaconda/anaconda3-2022b/
source ${CONDA_ROOT}/etc/profile.d/conda.sh
module load anaconda/2022b
conda activate terminator-nightly
echo OUTPUTDIR

python to_etab.py \
    --output_dir=OUTPUTDIR \
    --save_type=numpy \
    --num_cores=64