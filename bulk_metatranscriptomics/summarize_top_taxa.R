### a script for summarizing top phylogenetic levels for xtree, kraken, or metaphlan4 data

library(tidyverse)
library(ggtree)
library(umap)
library(ComplexHeatmap)
library(ggbeeswarm)
library(ComplexUpset)
library(tidytext)
library(circlize)
library(vegan)
library(ggraph)
library(broom)
library(reshape2)
library(igraph)
library(ggrepel)
library(broom)
library(ggpubr)
library(ape)
library(phytools)


### KRAKEN
brackout = read.csv('~/Dropbox (Mason Lab)/radx/data_packet_feb23/kraken_data_masked/wawa_bracken_merged_masked.tsv',sep='\t',header=T)
brackout_num = brackout %>% select(name,taxonomy_id,all_of(grep('num',colnames(.))))
colnames(brackout_num) = gsub('.masked.bracken_num','',colnames(brackout_num))
brackout_frac = brackout %>% select(name,taxonomy_id,all_of(grep('frac',colnames(.))))
colnames(brackout_frac) = gsub('.masked.bracken_frac','',colnames(brackout_frac))

toremove = brackout_num %>% filter(Negative.Control_Extraction_Batch06>quantile(brackout_num$Negative.Control_Extraction_Batch06,.75)) %>% select(taxonomy_id) %>% unlist %>% unname
brackout_num = brackout_num %>% filter(!(taxonomy_id %in% toremove))

toremove = brackout_frac %>% filter(Negative.Control_Extraction_Batch06>quantile(brackout_frac$Negative.Control_Extraction_Batch06,.75)) %>% select(taxonomy_id) %>% unlist %>% unname
brackout_frac = brackout_frac %>% filter(!(taxonomy_id %in% toremove))

# map tax id to 
taxmap = taxonomizr::getTaxonomy(brackout_num$taxonomy_id,sqlFile = '~/Dropbox (Mason Lab)/i4/i4_data_packet/accessionTaxa.sql') %>% data.frame %>%rownames_to_column('taxonomy_id') %>% filter(!(is.na(superkingdom))) %>% select(taxonomy_id,superkingdom,family,species) 
taxmap$taxonomy_id=as.numeric(taxmap$taxonomy_id)

brackout_frac = inner_join(brackout_frac,taxmap %>% select(taxonomy_id,superkingdom,family,species))%>% select(-taxonomy_id,-name) %>% melt %>% filter(variable!="Negative.Control_Extraction_Batch06",variable!='Positive.Control_Extraction_Batch06')
brackout_frac = brackout_frac %>% mutate(variable=gsub('_B06','',variable)) 

# bracken family
brackout_frac_bac = brackout_frac %>% filter(superkingdom == 'Bacteria') %>% select(-superkingdom) %>% group_by(family,variable) %>% summarise(value=sum(value))

# bracken vir family
brackout_frac_vir = brackout_frac %>% filter(superkingdom == 'Viruses') %>% select(-superkingdom,-family)

### METAPHLAN
metaph = read.table('~/Dropbox (Mason Lab)/radx/data_packet_feb23/metaphlan/merged_metaphlan_data.tsv',sep='\t',header=T) %>% filter(grepl('k__Bacteria',clade_name))%>% filter(grepl('f__',clade_name))%>% filter(!grepl('g__',clade_name)) %>% mutate(species = strsplit(clade_name,'f__') %>% map_chr(2) %>% gsub('_',' ',.)) %>% select(-clade_name)
colnames(metaph) = gsub('_metaphlan','',colnames(metaph)) 
toremove = metaph %>% filter(Negative.Control_Extraction_Batch06>quantile(metaph$Negative.Control_Extraction_Batch06,.75)) %>% select(species) %>% unlist %>% unname
metaph = metaph %>% filter(!(species %in% toremove))%>% melt %>% filter(variable!="Negative.Control_Extraction_Batch06",variable!='Positive.Control_Extraction_Batch06')

### XTREE BAC -- family
xtree_bac_family = read.delim('~/Dropbox (Mason Lab)/radx/data_packet_feb23/xtree/bacterial/GTDB_.01_.0.005_7_metatranscriptomics_ra.tsv')
xtree_bac_family = xtree_bac_family[rownames(xtree_bac_family)!='Unknown',] 
xtree_bac_family = xtree_bac_family[,colSums(xtree_bac_family)!=0]
xtree_bac_family = xtree_bac_family %>% t %>% data.frame(check.names=F) %>%  rownames_to_column('sample_id')
xtree_bac_family = xtree_bac_family %>% mutate(sample_id=gsub('_B06','',sample_id)) %>% column_to_rownames('sample_id') %>% t %>% data.frame %>% rownames_to_column('clade_name')
#xtree_bac_family = xtree_bac_family %>% mutate(species = strsplit(clade_name,';s__') %>% map_chr(2)) %>% select(-clade_name) 
toremove = xtree_bac_family %>% filter(Negative.Control_Extraction_Batch06>quantile(xtree_bac_family$Negative.Control_Extraction_Batch06,.75)) %>% select(clade_name) %>% unlist %>% unname
xtree_bac_family = xtree_bac_family %>% filter(!(clade_name %in% toremove))%>% melt %>% filter(variable!="Negative.Control_Extraction_Batch06",variable!='Positive.Control_Extraction_Batch06')
xtree_bac_family = xtree_bac_family %>% filter(grepl('f__',clade_name)) %>% mutate(clade_name = strsplit(clade_name,'f__') %>% map_chr(2) %>% strsplit(';g__') %>% map_chr(1))

### XTREE VIR
xtree_vir = read.delim('~/Dropbox (Mason Lab)/radx/data_packet_feb23/xtree/viral/vir_.2_.0.1_metatranscriptomics_species_ra.tsv')
xtree_vir = xtree_vir[rownames(xtree_vir)!='Unknown',] 
xtree_vir = xtree_vir[,colSums(xtree_vir)!=0]
xtree_vir = xtree_vir %>% t %>% data.frame(check.names=F) %>%  rownames_to_column('sample_id')
xtree_vir = xtree_vir %>% mutate(sample_id=gsub('_B06','',sample_id)) %>% column_to_rownames('sample_id') %>% t %>% data.frame %>% rownames_to_column('species')
#toremove = xtree_vir %>% filter(Negative.Control_Extraction_Batch06>quantile(xtree_vir$Negative.Control_Extraction_Batch06,.75)) %>% select(species) %>% unlist %>% unname
xtree_vir = xtree_vir %>% filter(!(species %in% toremove))%>% melt %>% filter(variable!="Negative.Control_Extraction_Batch06",variable!='Positive.Control_Extraction_Batch06')
#filter out partial genomes
xtree_vir = xtree_vir %>% filter(!grepl('segment',species),!grepl('hypothetical',species),!grepl('gene',species),!grepl('partial',species))

# load metadata
mdat = read.csv('~/Dropbox (Mason Lab)/radx/bulk/bulk_metadata.csv')# %>% filter(sample_category_wastewater == 'WEEKLY')
mdat$sampling.date = as.Date(mdat$sampling.date)
mdat$month = month.abb[as.numeric(format(mdat$sampling.date, "%m") )]
mdat$month = factor(mdat$month,levels = c('Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug'))
mdat$month2 = mdat$month
levels(mdat$month) = c(9,10,11,12,1,2,3,4,5,6,7,8)

locationmap = read.csv('~/Dropbox (Mason Lab)/radx/location_mapping_radx_TB.csv')
locationmap$NAME=toupper(locationmap$NAME)
mdat$site_name = toupper(mdat$site_name)
mdat = inner_join(mdat,locationmap,by=c('site_name'='NAME'))

mdat$month = as.numeric(as.character(mdat$month))

mdat = mdat %>% arrange(month)

mdat$month_norm = 2 * pi * mdat$month / max(mdat$month,na.rm=T)
mdat$month_cyclical = cos(mdat$month_norm)

mergedata <- function(data){
  data1 = data %>% filter(grepl('SL',variable))
  merged1 = inner_join(data1,mdat,by=c('variable'='ha_sample_id'))
  data2 = data %>% filter(!grepl('SL',variable)) %>% mutate(variable= gsub('\\.','-',variable))
  merged2 = inner_join(data2,mdat,by=c('variable'='sample_id'))
  merged1 = bind_rows(merged1,merged2) %>% filter(site_name_description!='n/a')
  return(merged1)
}


### metaph, xtree_bac_family, xtree_vir, brackout_frac_bac, brackout_frac_vir

alldata = bind_rows(metaph %>% dplyr::rename(taxa=species) %>% mutate(classifier = 'MetaPhlAn4'), xtree_bac_family%>% dplyr::rename(taxa=clade_name)%>% mutate(classifier = 'xtree/GTDB'), xtree_vir%>% dplyr::rename(taxa=species)%>% mutate(classifier = 'xtree/GenBank_viruses'), brackout_frac_bac%>% dplyr::rename(taxa=family)%>% mutate(classifier = 'kraken2/RefSeq_bacteria'), brackout_frac_vir%>% dplyr::rename(taxa=species)%>% mutate(classifier = 'kraken2/RefSeq_viruses'))

merged = mergedata(alldata)

#topfam = merged %>% group_by(LOCATION_TYPE,classifier,taxa) %>% summarise(mean = mean(value)) %>% slice_max(mean,n=25) %>% mutate(ordergroup = paste(LOCATION_TYPE,classifier))
#plot = ggplot(data = topfam %>% filter(!is.na(taxa))%>% filter(!grepl('vir',classifier)),aes(x = mean,y=taxa)) + geom_bar(stat='identity') + facet_grid(LOCATION_TYPE ~ classifier,scales='free') + theme(axis.text.x = element_text(angle=45,hjust=1)) 
#ggsave(plot = plot,'~/Dropbox (Mason Lab)/radx/plots/bacterial_family_summary_across_classifier.pdf',width=16,height=30)


wide = dcast(merged%>% filter(!grepl('vir',classifier)) %>% group_by(LOCATION_TYPE,classifier,taxa) %>% summarise(mean = mean(value)) ,taxa  ~  classifier,value.var = 'mean')
wide = wide %>% filter(!is.na(taxa)) %>% column_to_rownames('taxa')
wide[wide>0]=1
wide[wide!=1]=0
wide[is.na(wide)]= 0

pdf('~/Dropbox (Mason Lab)/radx/plots/bacterial_family_summary_across_classifier_upset.pdf',width=5,height=4)
upset(wide,sets=colnames(wide))
dev.off()


merged = merged %>% filter(grepl('vir',classifier)) %>% mutate(taxa =gsub(' complete sequence','',taxa))%>% mutate(taxa =gsub(' RNA 2','',taxa))%>% mutate(taxa =gsub(' RNA 3','',taxa))%>% mutate(taxa =gsub(', complete genome','',taxa))%>% mutate(taxa =gsub(' RNA 1','',taxa))%>% mutate(taxa =gsub(' RNA1','',taxa))%>% mutate(taxa =gsub(', genome','',taxa))


wide = dcast(merged%>% group_by(LOCATION_TYPE,classifier,taxa) %>% summarise(mean = mean(value)) ,taxa  ~  classifier,value.var = 'mean')
wide = wide %>% filter(!is.na(taxa)) %>% column_to_rownames('taxa')
wide[wide>0]=1
wide[wide!=1]=0
wide[is.na(wide)]= 0



pdf('~/Dropbox (Mason Lab)/radx/plots/viral_family_summary_across_classifier_upset.pdf',width=5,height=4)
upset(wide,sets=colnames(wide))
dev.off()

