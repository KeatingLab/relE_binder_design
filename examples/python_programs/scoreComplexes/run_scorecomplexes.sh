#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name scoreBinders
#SBATCH --tasks-per-node=1
#SBATCH --mem-per-cpu=32G #Request 4G of memory per CPU
#SBATCH -e scoreBinders.%j.err
#SBATCH -o scoreBinders.%j.out
#SBATCH -t 1-00:00:00
#SBATCH --gres=gpu:volta:1

export PYTHONPATH=$PYTHONPATH:/orcd/home/002/swans/relE_paper/relE_binder_design/
export PYTHONPATH=$PYTHONPATH:/orcd/home/002/swans/relE_paper/relE_binder_design/TERMinator_sscore

CUDA_LAUNCH_BLOCKING=1

dataset=/orcd/home/002/swans/relE_paper/relE_binder_design/examples/python_programs/scoreComplexes/7_connectFragments_relBhelix-2seeds_memfixed_topologyfixed_fusedDB.pdb
target_pdb=/orcd/home/002/swans/relE_paper/relE_binder_design/examples/python_programs/scoreComplexes/4FXE-ARG81-relax-noHyd_E_.pdb
TERMrepo=/orcd/home/002/swans/relE_paper/relE_binder_design/
model=/orcd/pool/008/swans/terminator_stuff/priority/230308_retrainCOORDinator_filtTrainData_etabnormpenalty_2
# model=/orcd/pool/008/swans/terminator_stuff/priority/230310_finetuneCOORDinator_sscore_mpnodeupdate
seqMode="consensus_aa"
scoreMode="interface_only"
custom_pep_reference=/orcd/pool/008/swans/proteome_amino_acid_distribution/iupred_humanproteome0.7cutoff_aaprob.pkl

SECONDS=0

poetry run python $TERMrepo/scripts/python/scoreBinders.py \
      --binder_dataset $dataset \
      --target_pdb $target_pdb \
      --model_dir $model \
      --seq_mode $seqMode \
      --score_mode $scoreMode \
      --custom_pep_reference $custom_pep_reference

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED
