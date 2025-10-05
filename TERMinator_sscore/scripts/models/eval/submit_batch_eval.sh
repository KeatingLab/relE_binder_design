#!/bin/bash
DTERMEN_ROCKLIN=/data1/groups/keatinglab/mlu/rocklin/rocklin_cleaned
FEATURES_ROCKLIN=/data1/groups/keatinglab/mlu/rocklin/rocklin_features
NRGTEN_ROCKLIN=/data1/groups/keatinglab/mlu/rocklin/rocklin_features
IN_ROCKLIN=/data1/groups/keatinglab/fosterb/rocklin_data/nrgten/rocklin.in
DTERMEN_MULTICHAIN=/data1/groups/keatinglab/fosterb/data/multichain_clean
FEATURES_MULTICHAIN=/data1/groups/keatinglab/fosterb/data/multichain_features
NRGTEN_MULTICHAIN=/data1/groups/keatinglab/fosterb/data/multichain_nrgten_small_motion_78
IN_MULTICHAIN=/data1/groups/keatinglab/fosterb/data/test_trim.in
DTERMEN_BCL2=/data1/groups/keatinglab/mlu/bcl2/clean_pdbs
FEATURES_BCL2=/data1/groups/keatinglab/mlu/bcl2/termless_features_with_sortcery
NRGTEN_BCL2=/data1/groups/keatinglab/fosterb/bcl2_data/nrgten
IN_BCL2=/data1/groups/keatinglab/fosterb/bcl2_data/nrgten/bcl2.in
DTERMEN_COVID=/data1/groups/keatinglab/fosterb/covid_data/clean_pdb
FEATURES_COVID=/data1/groups/keatinglab/fosterb/covid_data/features
NRGTEN_COVID=/data1/groups/keatinglab/fosterb/covid_data/nrgten
IN_COVID=/data1/groups/keatinglab/fosterb/covid_data/covid.in
for i in rocklin_results multichain_results bcl2_results covid_results
do
	for j in flex_ensembles_graph_cat_shared
	do
		MODEL_DIR=~/TERMinator_experiments/${j}
		OUTPUT_DIR=~/TERMinator_experiments/${j}/${i}
		if [ $i == bcl2_results ]
		then
			bash ~/TERMinator/scripts/models/eval/submit_eval.sh ${MODEL_DIR} ${FEATURES_BCL2} ${OUTPUT_DIR} ${IN_BCL2} ${NRGTEN_BCL2} ${DTERMEN_BCL2}
		elif [ $i == rocklin_results ]
		then
			bash ~/TERMinator/scripts/models/eval/submit_eval.sh ${MODEL_DIR} ${FEATURES_ROCKLIN} ${OUTPUT_DIR} ${IN_ROCKLIN} ${NRGTEN_ROCKLIN} ${DTERMEN_ROCKLIN}
		elif [ $i == covid_results ]
		then
			bash ~/TERMinator/scripts/models/eval/submit_eval.sh ${MODEL_DIR} ${FEATURES_COVID} ${OUTPUT_DIR} ${IN_COVID} ${NRGTEN_COVID} ${DTERMEN_COVID}
		else
			bash ~/TERMinator/scripts/models/eval/submit_eval.sh ${MODEL_DIR} ${FEATURES_MULTICHAIN} ${OUTPUT_DIR} ${IN_MULTICHAIN} ${NRGTEN_MULTICHAIN} ${DTERMEN_MULTICHAIN}
		fi
	done
done
# 0.18084923923015594 0.047488365322351456 0.005339636001735926 0.055779024958610535 0.05177825689315796 0.05326477810740471 0.15630945563316345 0.12791234254837036 0.1800616830587387 0.14251022040843964 0.19601035118103027 0.09344129264354706 0.09419145435094833 0.025428391993045807 0.11600253731012344

#/data1/groups/keatinglab/mlu/bcl2/termless_features_with_sortcery
#/data1/groups/keatinglab/mlu/rocklin/rocklin_features /home/gridsan/fbirnbaum/TERMinator_experiments/noise_batch_0/run_noise_level_${i}/rocklin_results /data1/groups/keatinglab/fosterb/rocklin_data/rocklin.in
# 	. /home/gridsan/fbirnbaum/TERMinator/scripts/models/eval/submit_eval.sh /home/gridsan/fbirnbaum/TERMinator_experiments/flex_ddp /data1/groups/keatinglab/mlu/bcl2/termless_features_sortcery_split /home/gridsan/fbirnbaum/TERMinator_experiments/flex_ddp/bcl2_results /data1/groups/keatinglab/mlu/bcl2/termless_features_sortcery_split/test.in
	#. /home/gridsan/fbirnbaum/TERMinator/scripts/models/eval/submit_eval.sh /home/gridsan/fbirnbaum/TERMinator_experiments/noise_batch_0/run_noise_level_${i} /data1/groups/keatinglab/fosterb/data/multichain_features /home/gridsan/fbirnbaum/TERMinator_experiments/noise_batch_0/run_noise_level_${i}/results /data1/groups/keatinglab/fosterb/data/multichain_features/test.in
 	#. /home/gridsan/fbirnbaum/TERMinator/scripts/models/eval/submit_eval.sh /home/gridsan/fbirnbaum/TERMinator_experiments/noise_batch_50/run_noise_level_0.5 /data1/groups/keatinglab/mlu/rocklin/rocklin_features /home/gridsan/fbirnbaum/TERMinator_experiments/noise_batch_50/run_noise_level_0.5/rocklin_results /data1/groups/keatinglab/fosterb/rocklin_data/rocklin.in
# noise_batch_0/run_noise_level_0.0 noise_batch_0/run_noise_level_0.05 noise_batch_0/run_noise_level_0.1 noise_batch_0/run_noise_level_0.15 noise_batch_0/run_noise_level_0.2 noise_batch_0/run_noise_level_0.5 noise_residue_batch_0/run_noise_level_0.0 noise_residue_batch_0/run_noise_level_0.05 noise_residue_batch_0/run_noise_level_0.1 noise_residue_batch_0/run_noise_level_0.15 noise_residue_batch_0/run_noise_level_0.2 noise_residue_batch_0/run_noise_level_0.5