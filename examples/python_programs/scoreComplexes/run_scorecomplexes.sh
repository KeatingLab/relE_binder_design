#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name scoreBinders
#SBATCH --tasks-per-node=1
#SBATCH --mem-per-cpu=32G #Request 4G of memory per CPU
#SBATCH -e scoreBinders.%j.err
#SBATCH -o scoreBinders.%j.out
#SBATCH -t 1-00:00:00
#SBATCH --gres=gpu:volta:1

source $(conda info --base)/etc/profile.d/conda.sh
conda activate terminator_sscore
export PYTHONPATH=$PYTHONPATH:~/local_code_mirror/peptide_binder_design
export PYTHONPATH=$PYTHONPATH:~/local_code_mirror/TERMinator_sscore

dataset=/data1/groups/keatinglab/swans/side_projects/lguan/sscore_propedia/1_lazyrename/posttrain_complexes_paths.txt
TERMrepo=/home/gridsan/sswanson/local_code_mirror/peptide_binder_design
model=/data1/groups/keatinglab/swans/TERMinator/230310_finetuneCOORDinator_sscore_mpnodeupdate/
seqMode="consensus_aa"
scoreMode="interface_only"
custom_pep_reference=/data1/groups/keatinglab/swans/proteome_amino_acid_distribution/iupred_humanproteome0.7cutoff_aaprob.pkl

SECONDS=0

python $TERMrepo/python/score_binders/scoreBinders.py \
      --complex_dataset $dataset \
      --model_dir $model \
      --seq_mode $seqMode \
      --score_mode $scoreMode \
      --custom_pep_reference $custom_pep_reference

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED
