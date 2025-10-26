#!/bin/bash
. /etc/profile.d/modules.sh

#######################################################################
# Set these variables before running
list_file=silentFromPDBs_denovoMCMC_list.txt
runfile=$PWD/run_relaxComplex.sh # the script that will be run for each element in the file
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
