#SBATCH --time=0-00:20:00

#log files:
#SBATCH -e run_multimer_jobs_%A_%a.err
#SBATCH -o run_multimer_jobs_%A_%a.out

#qos sets priority
# SBATCH --qos=low

# SBATCH -p gpu
#lower end GPUs might be sufficient for pairwise screens:
# SBATCH -C "gpu=2080Ti|gpu=3090"

#Reserve the entire GPU so no-one else slows you down
#SBATCH --gres=gpu:volta:1
#SBATCH --constraint=xeon-g6

#Limit the run to a single node
#SBATCH -N 1

#Adjust this depending on the node
#SBATCH --ntasks=8
#SBATCH --mem=64000

program=/home/gridsan/sswanson/local_code_mirror/AlphaPulldown/alphapulldown/run_multimer_jobs.py

echo "SLURM_ARRAY_TASK_ID: " $1
#mkdir -p output

SECONDS=0

python $program --mode=pulldown \
    --num_cycle=3 \
    --num_predictions_per_model=1 \
    --output_path=output \
    --data_dir=/data1/groups/keatinglab/alphafold_ss/dataset/ \
    --protein_lists=/data1/groups/keatinglab/swans/binderDesign_relE/analysis/round4/4FXE_relBhelixextensions_top5k-structscorebetterthanmedian-30ormoreres-seqstructscore_list.txt,relE.txt \
    --monomer_objects_dir=createFeatures/output/features,/data1/groups/keatinglab/swans/binderDesign_relE/data/alphafold_multimer/createIndividualFeatures_relE-ARG81/output/features \
    --no_pair_msa \
    --job_index $1

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED