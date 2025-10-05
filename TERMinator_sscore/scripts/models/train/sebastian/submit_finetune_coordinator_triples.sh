#!/bin/bash
. /etc/profile.d/modules.sh

# Set variables
TERMINATOR_REPO=/home/gridsan/sswanson/local_code_mirror/TERMinator_sscore

DATASET=/data1/groups/keatinglab/swans/TERMinator/finetune_multichain_features
MODEL_HPARAMS=$TERMINATOR_REPO/hparams/model/coordinator.json
RUN_HPARAMS=$TERMINATOR_REPO/hparams/run/finetune_sscore.json
RUNDIR=$PWD
TRAIN=train_filt_subsample.in
VALIDATION=validation_subsample.in
TEST=test_subsample.in
MASTER_PORT=1234
MPN_TYPE=dgat
NODE_RANK=$LLSUB_RANK

outDir=$PWD
export NUM_GPUS_PER_NODE=2
export NUM_NODES=$LLSUB_SIZE
export NUM_TRIALS=1

# Figure out whether this is the master process and if not, which node the master process is on
LLSUB_NAME="${LLSUB_OUTPUT_SUBDIR%/*}"
LLSUB_ID="${LLSUB_NAME#LLSUB.}"
FILE=master_node_${LLSUB_ID}.txt
echo $LLSUB_ID
MASTER_NODE=""
if [[ $LLSUB_RANK -eq 0 ]]; then
    # Main process
    MASTER_NODE=$SLURMD_NODENAME
    echo $MASTER_NODE > $FILE
    MASTER_ADDR=localhost
else
    # Worker process
    while [[ ! -e $FILE ]]; do
        echo "sleep 5"
        sleep 5
    done
    sleep 1 # just in case the file isn't ready yet
    MASTER_NODE=$(cat "$FILE")
    MASTER_ADDR=$MASTER_NODE
fi

# Record the variable values and call the training script
echo bash $TERMINATOR_REPO/scripts/models/train/run_train_gpu_ddp_triples.sh $NODE_RANK $MASTER_ADDR $DATASET $MODEL_HPARAMS $RUN_HPARAMS $RUNDIR $TERMINATOR_REPO
bash $TERMINATOR_REPO/scripts/models/train/run_train_gpu_ddp_triples.sh $NODE_RANK $MASTER_ADDR $DATASET $MODEL_HPARAMS $RUN_HPARAMS $RUNDIR $TERMINATOR_REPO