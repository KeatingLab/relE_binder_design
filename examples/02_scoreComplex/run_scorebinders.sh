#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name scoreBinders
#SBATCH --tasks-per-node=1
#SBATCH --mem-per-cpu=8G #Request 4G of memory per CPU
#SBATCH -e scoreBinders.%j.err
#SBATCH -o scoreBinders.%j.out
#SBATCH -t 1-00:00:00
#SBATCH --gres gpu:volta:1

export PYTHONPATH=$PYTHONPATH:/orcd/home/002/swans/relE_paper/relE_binder_design/
export PYTHONPATH=$PYTHONPATH:/orcd/home/002/swans/relE_paper/relE_binder_design/TERMinator_sscore

complex_dataset=2KC8_list.txt
target_chain_id="A"
binder_chain_id="B"
TERMrepo=/orcd/home/002/swans/relE_paper/relE_binder_design/
model=/orcd/pool/008/swans/terminator_stuff/priority/230308_retrainCOORDinator_filtTrainData_etabnormpenalty_2
seqMode="consensus_aa"
scoreMode="interface_only"
custom_pep_reference=/orcd/pool/008/swans/proteome_amino_acid_distribution/iupred_humanproteome0.7cutoff_aaprob.pkl

echo $PYTHONPATH

SECONDS=0

python $TERMrepo/scripts/python/scoreBinders.py \
      --complex_dataset $complex_dataset \
      --target_chain_id $target_chain_id \
      --binder_chain_id $binder_chain_id \
      --model_dir $model \
      --seq_mode $seqMode \
      --score_mode $scoreMode \
      --custom_pep_reference $custom_pep_reference

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED
