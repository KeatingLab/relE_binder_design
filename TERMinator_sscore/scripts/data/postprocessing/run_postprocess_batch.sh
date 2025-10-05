DTERMEN_DATA_ROOT_ROCKLIN=/data1/groups/keatinglab/mlu/rocklin/rocklin_cleaned
PDB_ROOT_ROCKLIN=/data1/groups/keatinglab/mlu/rocklin/rocklin_pdbs
DTERMEN_DATA_ROOT_MULTICHAIN=/data1/groups/keatinglab/fosterb/data/multichain_clean
PDB_ROOT_MULTICHAIN=/data1/groups/keatinglab/fosterb/data/multichain_raw
DTERMEN_DATA_ROOT_BCL2=/data1/groups/keatinglab/fosterb/bcl2_data/clean_pdbs
PDB_ROOT_BCL2=/data1/groups/keatinglab/mlu/bcl2/default_dtermen_actual_peptide
DTERMEN_DATA_ROOT_COVID=/data1/groups/keatinglab/fosterb/covid_data/clean_pdb
PDB_ROOT_COVID=/data1/groups/keatinglab/fosterb/covid_data/raw_pdb
# 0.05326477810740471 0.15630945563316345 0.12791234254837036 0.1800616830587387 0.14251022040843964 0.19601035118103027 0.09344129264354706 0.09419145435094833 0.025428391993045807 0.11600253731012344
# for i in 0.18084923923015594 0.047488365322351456 0.005339636001735926 0.055779024958610535 0.05177825689315796
SCRIPT=submit_etab_numpy.sh
for j in pairwise_dist noise_batch_1 noise_residue_batch_1 noise_torsion_batch_0.1_0.1 flex_ensembles_pos_variation
do
	for i in rocklin_results
	do
		OUTPUT_DIR=~/TERMinator_experiments/${j}/${i}
		if [ $i == bcl2_results ]
		then
			bash ~/TERMinator/scripts/data/postprocessing/${SCRIPT} ${DTERMEN_DATA_ROOT_BCL2} ${PDB_ROOT_BCL2} ${OUTPUT_DIR}
		elif [ $i == rocklin_results ]
		then
			bash ~/TERMinator/scripts/data/postprocessing/${SCRIPT} ${DTERMEN_DATA_ROOT_ROCKLIN} ${PDB_ROOT_ROCKLIN} ${OUTPUT_DIR}
		elif [ $i == covid_results ]
		then
			bash ~/TERMinator/scripts/data/postprocessing/${SCRIPT} ${DTERMEN_DATA_ROOT_COVID} ${PDB_ROOT_COVID} ${OUTPUT_DIR}
		else
			bash ~/TERMinator/scripts/data/postprocessing/${SCRIPT} ${DTERMEN_DATA_ROOT_MULTICHAIN} ${PDB_ROOT_MULTICHAIN} ${OUTPUT_DIR}
		fi
	done
done


# flex_ensembles_raw flex_ensembles_graph_stack_unique flex_ensembles_graph flex_ensembles_graph_stack_shared flex_ensembles_embedded

# noise_batch_0/run_noise_level_0.0 noise_batch_0/run_noise_level_0.05 noise_batch_0/run_noise_level_0.1 noise_batch_0/run_noise_level_0.15 noise_batch_0/run_noise_level_0.2 noise_batch_0/run_noise_level_0.5 noise_residue_batch_0/run_noise_level_0.0 noise_residue_batch_0/run_noise_level_0.05 noise_residue_batch_0/run_noise_level_0.1 noise_residue_batch_0/run_noise_level_0.15 noise_residue_batch_0/run_noise_level_0.2 noise_residue_batch_0/run_noise_level_0.5