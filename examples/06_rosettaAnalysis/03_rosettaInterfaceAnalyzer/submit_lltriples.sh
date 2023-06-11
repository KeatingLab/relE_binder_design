#!/bin/bash
. /etc/profile.d/modules.sh

#######################################################################
# Set these variables before running
list_file=1_relax_denovoMCMC-top5_list.txt
runfile=$PWD/run_interfaceAnalyzer.sh
#######################################################################

#LLSUB_SIZE=192
#LLSUB_RANK=191
echo $LLSUB_SIZE
echo $LLSUB_RANK

NUM=$(($LLSUB_RANK + 1))
silent_file=$(sed "${NUM}q;d" $list_file)
echo $silent_file

echo bash $runfile $silent_file
bash $runfile $silent_file
