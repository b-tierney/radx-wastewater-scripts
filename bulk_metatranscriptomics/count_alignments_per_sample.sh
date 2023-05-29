#!/bin/bash

# count number of alignments above 90 vs number of viruses per sample

while read p; do

name=$(echo $p | cut -f2 -d/ | sed 's/_rnaviral_assembled//g')
count=$(grep '>' $p | wc -l)
count2=$(cat "$p"_blastout | cut -f1 | sort | uniq | wc -l)
echo -e "$name""\t""$count""\t""$count2" >> viral_alignment_counts.tsv

done<soft_config