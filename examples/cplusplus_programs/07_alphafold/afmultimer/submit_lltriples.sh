#!/bin/bash

. /etc/profile.d/modules.sh
source $(conda info --base)/etc/profile.d/conda.sh

module load cuda/11.1
conda activate AlphaPulldown
export PYTHONPATH=~/local_code_mirror/AlphaPulldown:~/local_code_mirror/alphafold_custom

MAXRAM=$(echo `ulimit -m` '/ 1024.0'|bc)
GPUMEM=`nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits|tail -1`
export XLA_PYTHON_CLIENT_MEM_FRACTION=`echo "scale=3;$MAXRAM / $GPUMEM"|bc`
export TF_FORCE_UNIFIED_MEMORY='1'

RUNLIST_LEN=5000
batch_size=$(($(($RUNLIST_LEN/$LLSUB_SIZE))+1))
runfile=run_multimer.sh # the script that will be run for each element in the file

# get the batch boundaries
let start=$batch_size*$LLSUB_RANK
let next=$LLSUB_RANK+1
let next=$batch_size*$next
if [[ $next -gt $RUNLIST_LEN ]]
then
  let end=$RUNLIST_LEN
else
  let end=$next
fi

batch_size=

# run the batch
i=$start
while [[ $i -lt $end ]]
do
  echo job $i
  bash $runfile $i
  i=$(($i + 1))
done