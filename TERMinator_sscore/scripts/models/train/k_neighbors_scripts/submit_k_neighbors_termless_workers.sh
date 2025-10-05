#!/bin/bash
IFS=","
. /etc/profile.d/modules.sh
mkdir -p /home/gridsan/fbirnbaum/TERMinator_experiments/noise_ddp

export NUM_GPUS_PER_NODE=1
export NUM_NODES=2
export NUM_TRIALS=15
for args in \
	dgat,1,d-13-11-1; \
	do set -- $args
	echo $1
	. /home/gridsan/fbirnbaum/TERMinator/scripts/models/train/submit_train_ddp.sh /data1/groups/keatinglab/fosterb/data/multichain_features /home/gridsan/fbirnbaum/TERMinator/hparams/model/coordinator.json /home/gridsan/fbirnbaum/TERMinator/hparams/run/k_neighbors_optuna.json /home/gridsan/fbirnbaum/TERMinator_experiments/k_neighbors_ddp /home/gridsan/fbirnbaum/TERMinator_experiments/k_neighbors_ddp /home/gridsan/fbirnbaum/TERMinator_experiments/k_neighbors_ddp 240 train validation test $1 $2 $3
done