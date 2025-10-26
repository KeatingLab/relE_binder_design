#!/bin/bash
  
#########################################
# Set these values
nameandsequencecsv=final_library_dnaseqs.txt
run_tisigner=run_tisigner.sh

# After the job is completed, you will need to parse the output to extract the DNA sequences
#########################################

#LLSUB_SIZE=1
#LLSUB_RANK=0

# Skip through the file, taking steps equal to the size of the batch
n_lines=$(grep -c ^ "$nameandsequencecsv")
step_size=$LLSUB_SIZE

echo "N lines: " $n_lines
echo "step size: "$step_size

# sed lines are 1-indexed
i=$(( $LLSUB_RANK + 1))
while [[ $i -le $n_lines ]]; do
    echo job $i
        
    line=$(sed "${i}q;d" $nameandsequencecsv)

    name=$(echo $line | cut -d "," -f 1)
    aasequence=$(echo $line | cut -d "," -f 2)

    echo "Name: "$name", Sequence: "$aasequence

    $run_tisigner $aasequence $name

    # increment
    i=$(($i + $step_size))
done

echo "Done!"