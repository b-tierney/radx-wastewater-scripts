#! /bin/bash -l
#SBATCH --job-name=viralmetaspades
#SBATCH --output=logs/viralmetaspades_%A_%a.log
#SBATCH --time=168:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --partition=panda

conda activate viral_spades

samplefile=$1
sample=$(cat  $samplefile | cut -f1 |sed -n "${SLURM_ARRAY_TASK_ID}p")
read1_path=$(cat  $samplefile | cut -f2 |sed -n "${SLURM_ARRAY_TASK_ID}p")
read2_path=$(cat  $samplefile | cut -f3 |sed -n "${SLURM_ARRAY_TASK_ID}p")

rnaspades.py  -1 $read1_path -2 $read2_path -o "$sample"_rnaviral_assembled -t 24 -m 164
