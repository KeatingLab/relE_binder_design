#!/bin/bash
IFS=","
. /etc/profile.d/modules.sh
mkdir -p /home/gridsan/fbirnbaum/TERMinator_experiments/noise_batch_02
export NUM_GPUS_PER_NODE=1
export NUM_NODES=3
export NUM_TRIALS=1
for args in \
	dgat,0,localhost; \
	do set -- $args
	. /home/gridsan/fbirnbaum/TERMinator/scripts/models/train/submit_train_ddp.sh /data1/groups/keatinglab/fosterb/data/multichain_features /home/gridsan/fbirnbaum/TERMinator/hparams/model/noise_batch.json /home/gridsan/fbirnbaum/TERMinator/hparams/run/noise_optuna.json /home/gridsan/fbirnbaum/TERMinator_experiments/noise_batch_02 /home/gridsan/fbirnbaum/TERMinator_experiments/noise_batch_02 /home/gridsan/fbirnbaum/TERMinator_experiments/noise_batch_02 300 train validation test $1 $2 $3
done
