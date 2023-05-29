#! /bin/bash -l
#SBATCH --job-name=xtreeviral
#SBATCH --output=logs/xtreeviral-%A_%a.log
#SBATCH --time=168:00:00
#SBATCH --mem=224G
#SBATCH --cpus-per-task=12
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

mkdir -p prof
mkdir -p prof_bac

#cat $read1 $read2 |  ~/UTree/xtree_recompile --seqs - --threads 12 --db /athena/masonlab/scratch/users/btt4001/databases/xtree_db/xtree_db_genbankviralcomplete_99percclust_k17_comp0_20220726/genbankviralcomplete_99percclust_k17_comp0_20220726.xtr --ref-out prof/"$samplename".ref --cov-out prof/"$samplename".cov --redistribute


cat $read1 $read2 |  ~/UTree/xtree_recompile --seqs - --threads 12 --db /athena/masonlab/scratch/users/btt4001/databases/xtree_db/xtree_db_gtdb207_kmer29_comp2_20220722/GTDB207_filt.xtr --ref-out prof_bac/"$samplename".ref --cov-out prof_bac/"$samplename".cov --redistribute


cat $read1 $read2 |  ~/UTree/xtree_recompile --seqs - --threads 12 --db /athena/masonlab/scratch/users/btt4001/databases/xtree_db/xtree_db_refseq_viral_kmer17_comp0_20220726/refseq_viral_kmer17_comp0_20220726.xtr --ref-out prof/"$samplename".ref --cov-out prof/"$samplename".cov --redistribute

cat $read1 $read2 |  ~/UTree/xtree_recompile --seqs - --threads 12 --db /athena/masonlab/scratch/users/btt4001/databases/xtree_db/xtree_db_genbank_fungi_protozoa_complete_20220726/genbank_fungi_protozoa_complete_20220726.xtr --ref-out prof_fp/"$samplename".ref --cov-out prof_fp/"$samplename".cov --redistribute


