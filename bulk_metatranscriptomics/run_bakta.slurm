#! /bin/bash -l
#SBATCH --job-name=bakta
#SBATCH --output=logs/bakta_%A_%a.log
#SBATCH --time=168:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=24
#SBATCH --partition=panda_physbio

conda activate bakta

samplefile=$1
basename=$(cat  $samplefile | cut -f1 |sed -n "${SLURM_ARRAY_TASK_ID}p")
seqfile=$(cat  $samplefile | cut -f2 |sed -n "${SLURM_ARRAY_TASK_ID}p")

mkdir -p "$basename"_tmp

bakta --threads 24 --tmp-dir "$basename"_tmp  --skip-pseudo --prefix "$basename"_bakta --db /athena/masonlab/scratch/users/btt4001/databases/bakta_db/db $seqfile
