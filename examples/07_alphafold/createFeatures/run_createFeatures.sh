#!/bin/bash

#script copied from README.MD https://github.com/KosinskiLab/AlphaPulldown/blob/main/example_1.md

#A typical run takes couple of hours but may be much longer
#SBATCH --job-name=alphafold_array
#SBATCH --time=1:00:00

#log files:
#SBATCH -e create_individual_features_%A_%a.err
#SBATCH -o create_individual_features_%A_%a.out

#qos sets priority
# SBATCH --qos=low

#Limit the run to a single node
#SBATCH -N 1

#Adjust this depending on the node
#SBATCH --ntasks=8
#SBATCH --mem=32G

#module load HMMER/3.3.2-gompic-2020b
#module load HH-suite/3.3.0-gompic-2020b
#module load Anaconda3
source $(conda info --base)/etc/profile.d/conda.sh
conda activate AlphaPulldown
export PYTHONPATH=~/local_code_mirror/AlphaPulldown:~/local_code_mirror/alphafold_custom
program=/home/gridsan/sswanson/local_code_mirror/AlphaPulldown/alphapulldown/create_individual_features.py

mkdir -p output/features

SECONDS=0

python $program \
  --fasta_paths=/data1/groups/keatinglab/swans/binderDesign_relE/analysis/round4/4FXE_relBhelixextensions_top5k-structscorebetterthanmedian-30ormoreres-seqstructscore.fasta \
  --data_dir=/data1/groups/keatinglab/alphafold_ss/dataset/ \
  --save_msa_files=True \
  --output_dir=${PWD}/output/features \
  --use_precomputed_msas=False \
  --max_template_date=2050-01-01 \
  --skip_existing=True \
  --no_msa_or_templates=True
  #  --seq_index=$SLURM_ARRAY_TASK_ID

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED