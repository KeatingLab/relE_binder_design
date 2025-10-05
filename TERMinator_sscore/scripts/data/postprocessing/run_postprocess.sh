DTERMEN_DATA_ROOT=~/keatinglab_shared/fosterb/data/multichain_clean
PDB_ROOT=~/keatinglab_shared/fosterb/data/multichain_raw
OUTPUT_DIR=~/TERMinator_experiments/noise_ddp/run_noise_0.047488365322351456/results

bash ~/TERMinator/scripts/data/postprocessing/submit_etab.sh ${DTERMEN_DATA_ROOT} ${PDB_ROOT} ${OUTPUT_DIR}