#! /bin/bash -l
#SBATCH --job-name=blastvir
#SBATCH --output=logs/blastvir_%A_%a.log
#SBATCH --time=168:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=panda_physbio
conda activate
samplefile=$1
sample=$(cat  $samplefile | cut -f1 |sed -n "${SLURM_ARRAY_TASK_ID}p")

blastn -query $sample -db refseq_viral_blastdb -out "$sample"_blastout -outfmt "6 std qlen slen qcovs" -max_target_seqs 25000 -perc_identity 90 -num_threads 4
