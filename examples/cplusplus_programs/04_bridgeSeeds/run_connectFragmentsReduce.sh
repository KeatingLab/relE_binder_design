#!/bin/bash

bin=/home/gridsan/sswanson/local_code_mirror/peptide_binder_design/cplusplus/bin/connectFragmentsReduce

SECONDS=0

echo $bin 
$bin

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED