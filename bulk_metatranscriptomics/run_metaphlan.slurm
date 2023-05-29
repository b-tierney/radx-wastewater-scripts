#! /bin/bash -l
#SBATCH --job-name=xtreeviral
#SBATCH --output=logs/metaphlan-%A_%a.log
#SBATCH --time=168:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=panda_physbio

### set up variables and import data
samplefile=$1

samplefile=$1
samplename=$(cat $samplefile | cut -f1 |sed -n "${SLURM_ARRAY_TASK_ID}p")
read1=$(cat  $samplefile | cut -f2 |sed -n "${SLURM_ARRAY_TASK_ID}p")
read2=$(cat  $samplefile | cut -f3 |sed -n "${SLURM_ARRAY_TASK_ID}p")
echo $samplename
echo $read1
echo $read2
conda activate metaphlan4
metaphlan "${read1}","${read2}" --input_type fastq -o "${samplename}"_metaphlan --nproc 8 --bowtie2db /athena/masonlab/scratch/users/btt4001/databases/metaphlan4_20230119/ --bowtie2out "${samplename}"_bowtie2out


