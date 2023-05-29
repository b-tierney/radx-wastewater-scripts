# merge all checkv/blast quality data and add sample names (for wastewater)

files = read.csv('quality_data_checkv_locs',header=F) %>% unlist %>% unname

data=list()
for(f in files){
print(f)
	data[[f]] = read.table(f,sep='\t',header=T) %>% mutate(Sample = strsplit(f,'/') %>% map_chr(2))
}

data_merged1 = bind_rows(data)

write.csv(data_merged1,'quality_summary_checkv_all.tsv',quote=F)


files = read.csv('blastout_locs',header=F) %>% unlist %>% unname

data=list()
for(f in files){
print(f)
	data[[f]] = read.table(f,sep='\t',header=F) %>% mutate(Sample = strsplit(f,'/') %>% map_chr(2))
}

data_merged2 = bind_rows(data) %>% dplyr::rename(contig_id=V1)

write.csv(data_merged2,'quality_summary_blast_all.tsv',quote=F)

out = inner_join(data_merged1,data_merged2,by=c('contig_id','Sample'))

saveRDS(out,'quality_alignment_data_aligned.rds')

mapping=read.delim('refseq_name_mapping',header=F)

merged = left_join(out,mapping,by=c('V2'='V1'))
merged = merged %>% dplyr::rename(viral_name = V2.y)

write.table(merged,'quality_alignment_data_aligned_viral_data.tsv',sep='\t',quote=F)






