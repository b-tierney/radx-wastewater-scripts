#! /bin/bash -l
#SBATCH --job-name=qcAndKraken
#SBATCH --output=logs/qcAndKraken-%A_%a.log
#SBATCH --time=168:00:00
#SBATCH --mem=96G
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
#/home/btt4001/kraken2-2.1.2/kraken2 --db /home/cem2009/ReferenceFiles/kraken2/databases/humanmicroorg --threads 8 --report ${samplename}.raw.k2.report --output ${samplename}.raw.k2.output --paired $read1 $read2
#/home/btt4001/Bracken-2.6.2/bracken -r 150 -d  /home/cem2009/ReferenceFiles/kraken2/databases/humanmicroorg -l S -i $samplename.raw.k2.report -o $samplename.raw.bracken


### clumpify (rapidly group overlapping reads into clumps)
echo '///////////////////////////////////////////////// STARTING QC' >> "$samplename"_log
clumpify.sh in=$read1 \
            in2=$read2 \
            -Xmx60G \
            out="${samplename}"_R1_clumped.fastq.gz \
            out2="${samplename}"_R2_clumped.fastq.gz \
            dedupe=t \
            dupesubs=2 \
            qin=33 \
            threads=8 \
            overwrite=true \
            optical=f 2>> "$samplename"_log
echo '///////////////////////////////////////////////// CLUMPIFY DONE' >> "$samplename"_log
echo; echo


### bbduk (deduplicate using kmers)
echo '///////////////////////////////////////////////// STARTING BBDUK' >> "$samplename"_log
bbduk.sh in="${samplename}"_R1_clumped.fastq.gz \
           -Xmx60G \
           in2="${samplename}"_R2_clumped.fastq.gz  \
           out="${samplename}"_R1_clumped_bbduk.fastq.gz \
           out2="${samplename}"_R2_clumped_bbduk.fastq.gz \
           interleaved=f \
           stats="${samplename}"_stats_bbduk \
           overwrite=true \
           qout=33 \
           trd=t \
           hdist=1 \
           k=27 \
           ktrim="r" \
           mink=8 \
           overwrite=true \
           trimq=10 \
           qtrim='rl' \
           threads=8 \
           minlength=51 \
           maxns=-1 \
           minbasefrequency=0.05 \
           ecco=f \
           prealloc=t \
           ref=/athena/masonlab/scratch/users/btt4001/databases/atlas_databases/adapters.fa 2>> "$samplename"_log
rm -rf "${samplename}"_R1_clumped.fastq.gz "${samplename}"_R2_clumped.fastq.gz
echo '///////////////////////////////////////////////// BBDUK DONE' >> "$samplename"_log
echo; echo


### bowtie on human reads
echo '///////////////////////////////////////////////// STARTING BOWTIE' >> "$samplename"_log
bowtie2 -p 8 \
        -x /athena/masonlab/scratch/users/btt4001/databases/hg38/GRCh38_noalt_as/GRCh38_noalt_as \
        -1 "${samplename}"_R1_clumped_bbduk.fastq.gz \
        -2 "${samplename}"_R2_clumped_bbduk.fastq.gz \
        --al-conc-gz "${samplename}"_human_reads \
        --very-sensitive-local \
        --un-conc-gz "${samplename}"_SAMPLE_hg38removed > "${samplename}"_SAMPLE_mapped_and_unmapped.sam 2>> "$samplename"_log

mv "${samplename}"_SAMPLE_hg38removed.1 "${samplename}"_R1_clumped_bbduk_hg38removed.fastq.gz
mv "${samplename}"_SAMPLE_hg38removed.2 "${samplename}"_R2_clumped_bbduk_hg38removed.fastq.gz

sambamba view -t 8 -S "${samplename}"_SAMPLE_mapped_and_unmapped.sam -f bam -o "$samplename".hg38.bam
rm -rf "${samplename}"_R1_clumped_bbduk.fastq.gz "${samplename}"_R2_clumped_bbduk.fastq.gz
rm -rf "${samplename}"_SAMPLE_mapped_and_unmapped.sam
echo '///////////////////////////////////////////////// BOWTIE DONE' >> "$samplename"_log
echo; echo


### repair remaining reads
echo '///////////////////////////////////////////////// STARTING REPAIRING' >> "$samplename"_log
repair.sh in="${samplename}"_R1_clumped_bbduk_hg38removed.fastq.gz \
          in2="${samplename}"_R2_clumped_bbduk_hg38removed.fastq.gz \
          out1="${samplename}"_R1_clumped_bbduk_hg38removed_repaired.fastq.gz \
          out2="${samplename}"_R2_clumped_bbduk_hg38removed_repaired.fastq.gz 2>> "$samplename"_log

rm -rf "${samplename}"_R1_clumped_bbduk_hg38removed.fastq.gz "${samplename}"_R2_clumped_bbduk_hg38removed.fastq.gz
echo '///////////////////////////////////////////////// REPAIRING DONE' >> "$samplename"_log
echo; echo


### tadpole (kmer-based assembler)
echo '///////////////////////////////////////////////// STARTING TADPOLE' >> "$samplename"_log
tadpole.sh prealloc=1 \
           -Xmx60G \
           in="${samplename}"_R1_clumped_bbduk_hg38removed_repaired.fastq.gz \
           in2="${samplename}"_R2_clumped_bbduk_hg38removed_repaired.fastq.gz \
           out="${samplename}"_R1_clumped_bbduk_hg38removed_repaired_tadpole.fastq.gz \
           out2="${samplename}"_R2_clumped_bbduk_hg38removed_repaired_tadpole.fastq.gz \
           mode=correct \
           threads=8 \
           ecc=t \
           ecco=t 2>> "$samplename"_log

rm -rf "${samplename}"_R1_clumped_bbduk_hg38removed_repaired.fastq.gz "${samplename}"_R2_clumped_bbduk_hg38removed_repaired.fastq.gz
echo '///////////////////////////////////////////////// TADPOLE DONE' >> "$samplename"_log


### kraken/bracken on final data
/home/btt4001/kraken2-2.1.2/kraken2 --confidence .2  --db /home/cem2009/ReferenceFiles/kraken2/databases/humanmicroorg --threads 8 --report ${samplename}.qc.k2.report --output ${samplename}.qc.k2.output --paired "${samplename}"_R1_clumped_bbduk_hg38removed_repaired_tadpole.fastq.gz "${samplename}"_R2_clumped_bbduk_hg38removed_repaired_tadpole.fastq.gz
/home/btt4001/Bracken-2.6.2/bracken -r 150  -d   /home/cem2009/ReferenceFiles/kraken2/databases/humanmicroorg -r 150 -l S -i $samplename.qc.k2.report -o $samplename.qc.bracken


### export
echo [`date`] Finished job

