#!/bin/bash
  
#########################################
# Set these values
nameandsequencecsv=final_library_aaseqs.txt
template_file=DNAWORKS_template.inp
dnaworks=DNAWorks/dnaworks
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

    # Create a file where each line contains:
    # NAME AASEQUENCE
    # note: NAME must be unique

    # Read file, go to specific line
    # Split line by , and get both values

    line=$(sed "${i}q;d" $nameandsequencecsv)

    name=$(echo $line | cut -d "," -f 1)
    aasequence=$(echo $line | cut -d "," -f 2)

    echo "Name: "$name", Sequence: "$aasequence

    # Copy the template file, replace the values with sed
    new_input_file=${name}.input
    if [ -f $new_input_file ]; then
        echo "File already exists: job was already run, or NAME values are not unique"
    else
        cp $template_file $new_input_file

        sed -i -e "s/NAME/$name/g" $new_input_file
        sed -i -e "s/AASEQUENCE/$aasequence/g" $new_input_file

        # Call DNAWorks
        $dnaworks $new_input_file

    fi

    # increment
    i=$(($i + $step_size))
done

echo "Done!"