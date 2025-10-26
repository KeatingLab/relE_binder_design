#!/bin/bash

bin=../../fragment_tools_cpp/connectFragmentsReduce

SECONDS=0

echo $bin 
$bin

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED