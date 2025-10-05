#!/bin/bash
IFS=","
. /etc/profile.d/modules.sh
mkdir -p /home/gridsan/fbirnbaum/TERMinator_experiments/flex_ensembles_test

export NUM_GPUS_PER_NODE=1
export NUM_NODES=2
export NUM_TRIALS=1
for args in \
	dgat,1,d-14-6-1; \
	do set -- $args
	echo $1
	. /home/gridsan/fbirnbaum/TERMinator/scripts/models/train/submit_train_ddp.sh /data1/groups/keatinglab/fosterb/data/multichain_features /home/gridsan/fbirnbaum/TERMinator/hparams/model/flex.json /home/gridsan/fbirnbaum/TERMinator/hparams/run/flex.json /home/gridsan/fbirnbaum/TERMinator_experiments/flex_ensembles_test /home/gridsan/fbirnbaum/TERMinator_experiments/flex_ensembles_test /home/gridsan/fbirnbaum/TERMinator_experiments/flex_ensembles_test 240 train validation test $1 $2 $3
done