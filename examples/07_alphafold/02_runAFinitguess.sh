#!/bin/bash

#SBATCH --job-name=AFinitialguesss
#SBATCH --time=1-00:00:00

#log files:
#SBATCH -e run_AFinitialguesss_%A.err
#SBATCH -o run_AFinitialguesss_%A.out

#Reserve the entire GPU so no-one else slows you down
#SBATCH --gres=gpu:volta:1
#SBATCH --constraint=xeon-g6

#Adjust this depending on the node
#SBATCH --ntasks=8
#SBATCH --mem=180G

module load cuda/11.1
conda activate dl_binder_design

#MAXRAM=$(echo `ulimit -m` '/ 1024.0'|bc)
#GPUMEM=`nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits|tail -1`
#export XLA_PYTHON_CLIENT_MEM_FRACTION=`echo "scale=3;$MAXRAM / $GPUMEM"|bc`
#export TF_FORCE_UNIFIED_MEMORY='1'

dl_binder_design_path=dl_binder_design
silent=$PWD/my_designs.silent
af_dir=alphafold/dataset/

echo $LLSUB_SIZE
echo $LLSUB_RANK

SECONDS=0

$dl_binder_design_path/af2_initial_guess/interfaceAF2predict.py -silent $silent \
	-af_dir $af_dir

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED
