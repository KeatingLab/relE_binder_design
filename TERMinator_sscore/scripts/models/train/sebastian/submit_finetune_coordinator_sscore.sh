#!/bin/bash


IFS=","
. /etc/profile.d/modules.sh
TERMinatorRepo=/home/gridsan/sswanson/local_code_mirror/TERMinator_sscore
outRepo=$PWD
mkdir -p $outRepo
export NUM_GPUS_PER_NODE=1
export NUM_NODES=$LLSUB_SIZE
export NUM_TRIALS=1
for args in \
        dgat,0,localhost; \
        do set -- $args
        . $TERMinatorRepo/scripts/models/train/submit_train_ddp.sh /data1/groups/keatinglab/swans/TERMinator/finetune_multichain_features $TERMinatorRepo/hparams/model/coordinator.json  $TERMinatorRepo/hparams/run/finetune_sscore.json $outRepo $outRepo $outRepo 24 train_filt_subsample validation_subsample test_subsample $1 $2 $3
done