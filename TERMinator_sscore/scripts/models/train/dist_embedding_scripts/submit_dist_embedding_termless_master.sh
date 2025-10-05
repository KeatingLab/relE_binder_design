#!/bin/bash
IFS=","
. /etc/profile.d/modules.sh
mkdir -p /data1/groups/keatinglab/fosterb/COORDinator_experiments/pairwise_dist_filt_ddp
export NUM_GPUS_PER_NODE=1
export NUM_NODES=3
export NUM_TRIALS=1
for args in \
	dgat,0,localhost; \
	do set -- $args
	. /home/gridsan/fbirnbaum/TERMinator/scripts/models/train/submit_train_ddp.sh /data1/groups/keatinglab/fosterb/data/multichain_features /home/gridsan/fbirnbaum/TERMinator/hparams/model/pairwise_dist.json /home/gridsan/fbirnbaum/TERMinator/hparams/run/pairwise_dist.json /data1/groups/keatinglab/fosterb/COORDinator_experiments/pairwise_dist_filt_ddp /data1/groups/keatinglab/fosterb/COORDinator_experiments/pairwise_dist_filt_ddp /data1/groups/keatinglab/fosterb/COORDinator_experiments/pairwise_dist_filt_ddp 24 train_filt validation test $1 $2 $3
done
