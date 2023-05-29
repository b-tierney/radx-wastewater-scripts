#! /bin/bash -l
#SBATCH --job-name=checkv
#SBATCH --output=logs/checkv_%A_%a.log
#SBATCH --time=168:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
#SBATCH --partition=panda_physbio

conda activate viral_spades

samplefile=$1
sample=$(cat  $samplefile | cut -f1 |sed -n "${SLURM_ARRAY_TASK_ID}p")
echo $sample
checkv end_to_end $sample "$sample"_checkv -t 4 -d /athena/masonlab/scratch/users/btt4001/checkv_database/checkv-db-v1.0
