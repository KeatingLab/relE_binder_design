#!/bin/bash
#SBATCH -N 1
#SBATCH --time=30:00:00
#SBATCH -o gen_MSA_QUERY.out
#SBATCH -e gen_MSA_QUERY.err


hhblits -i query.seq -d /data1/groups/keatinglab/fosterb/sequence_data/uniclust30/uniclust30_2018_08 -oa3m QUERY.a3m -cpu 32 -n 1