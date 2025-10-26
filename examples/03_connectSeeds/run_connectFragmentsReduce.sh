#!/bin/bash

# This tool is only necessary if you use run the mapper in parallel and want to merge files

bin=../../fragment_tools_cpp/connectFragmentsReduce

SECONDS=0

echo $bin 
$bin

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED
