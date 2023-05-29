#! /bin/bash -l
#SBATCH --job-name=qcAndKraken
#SBATCH --output=logs/qcAndKraken-%A_%a.log
#SBATCH --time=168:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=panda_physbio

# qc paired end, non interleaved metagenomes with bbtools
# you'll arguments are the base file name (eg sample_1) followed by the two paired end files.
# note absolute paths in bowtie2 and bbduk commands.
# adapters.fa file comes from the "atlas metagenome" workflow, available on github
# by default takes 10 threads + ~125G of RAM

# ----------------------------------------------------------------------------------------
# changelog
#
# 2022-01-19:
# incorporated kraken/bracken steps to run before/after QC steps
# changed $read2_path definition for HudsonAlpha naming
# changed $sample definition to strip barcode info (can always re-connect if needed)
# saving hg38 reads (converting sam to bam)
#
# ----------------------------------------------------------------------------------------

### initiate
echo [`date`] Started job
conda activate /home/btt4001/miniconda3/envs/atlasenv/

### set up variables and import data
samplefile=$1

samplefile=$1
samplename=$(cat $samplefile | cut -f1 |sed -n "${SLURM_ARRAY_TASK_ID}p")
read1=$(cat  $samplefile | cut -f2 |sed -n "${SLURM_ARRAY_TASK_ID}p")
read2=$(cat  $samplefile | cut -f3 |sed -n "${SLURM_ARRAY_TASK_ID}p")
echo $samplename

### kraken/bracken on raw data
/home/btt4001/kraken2-2.1.2/kraken2 --confidence .2 --db /home/cem2009/ReferenceFiles/kraken2/databases/humanmicroorg --threads 8 --report ${samplename}.masked.k2.report --output ${samplename}.masked.k2.output --paired $read1 $read2
/home/btt4001/Bracken-2.6.2/bracken -r 150 -d  /home/cem2009/ReferenceFiles/kraken2/databases/humanmicroorg -l S -i $samplename.masked.k2.report -o $samplename.masked.bracken



