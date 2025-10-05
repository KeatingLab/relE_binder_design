#!/bin/bash
#SBATCH -N 1
#SBATCH --time=30:00:00
#SBATCH -o design_ID.out
#SBATCH -e design_ID.err

bin=/data1/groups/keatinglab/MST_workspace/MST/bin/enerTable
cd PDBDIR

/data1/groups/keatinglab/MST_workspace/MST/bin/enerTable --e ETAB \
	--opt --lc 1 --s SEQUENCE
	
