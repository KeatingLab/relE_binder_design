#!/bin/bash
IFS=","
. /etc/profile.d/modules.sh
mkdir -p /home/gridsan/fbirnbaum/TERMinator_experiments/noise_residue_batch_50

export NUM_GPUS_PER_NODE=1
export NUM_NODES=3
export NUM_TRIALS=1
for args in \
	dgat,1,d-14-16-1 \
	dgat,2,d-14-16-1; \
	do set -- $args
	echo $1
	. /home/gridsan/fbirnbaum/TERMinator/scripts/models/train/submit_train_ddp.sh /data1/groups/keatinglab/fosterb/data/multichain_features /home/gridsan/fbirnbaum/TERMinator/hparams/model/noise_residue_batch.json /home/gridsan/fbirnbaum/TERMinator/hparams/run/noise_optuna.json /home/gridsan/fbirnbaum/TERMinator_experiments/noise_residue_batch_50 /home/gridsan/fbirnbaum/TERMinator_experiments/noise_residue_batch_50 /home/gridsan/fbirnbaum/TERMinator_experiments/noise_residue_batch_50 240 train validation test $1 $2 $3
done