

library(tidyverse)
library(umap)
library(ComplexHeatmap)
library(circlize)
library(vegan)
library(ggraph)
library(broom)
library(reshape2)
library(ggpubr)
library(igraph)
library(qdapTools)
library(ggrepel)
library(ggbeeswarm)

setwd('~/Dropbox (Mason Lab)/radx/data_packet_feb23/abricate/')

readcountspostqc = read.delim('../raw_sample_readcounts.tsv',header=F) %>% mutate(sample_id=gsub('_B06','',V1))%>% mutate(V1= gsub('\\.','-',V1)) %>% select(-V1)

data = read.delim('abricate_vfdb.tab',check.names=F) %>% mutate(sample_id = strsplit(`#FILE`,'/') %>% map_chr(2) %>% gsub('_rnaviral_assembled','',.))
data = data %>% mutate(sample_id=gsub('_B06','',sample_id))%>% mutate(sample_id= gsub('\\.','-',sample_id))
#data_wide = data %>% dcast(sample_id ~ GENE) 

datatemp = data %>% dcast(GENE + RESISTANCE ~ sample_id) 
res=mtabulate(strsplit(as.character(datatemp$RESISTANCE), ";")) 
res = bind_cols(datatemp %>% select(GENE),res) 
res = inner_join(datatemp %>% select(-RESISTANCE) %>% melt %>% filter(value ==1) %>% select(-value),res)

# load metadata

mdat = read.csv('../../bulk/bulk_metadata.csv')# %>% filter(sample_category_wastewater == 'WEEKLY')
mdat$sampling.date = as.Date(mdat$sampling.date)
mdat$month = month.abb[as.numeric(format(mdat$sampling.date, "%m") )]
mdat$month = factor(mdat$month,levels = c('Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug'))
mdat$month2 = mdat$month
levels(mdat$month) = c(9,10,11,12,1,2,3,4,5,6,7,8)

locationmap = read.csv('../../location_mapping_radx_TB.csv')
locationmap$NAME = toupper(locationmap$NAME)
mdat$site_name = toupper(mdat$site_name)
mdat = left_join(mdat,locationmap,by=c('site_name'='NAME'))
mdat$sample_description = gsub(' ','',mdat$sample_description)
mdat$month = as.numeric(as.character(mdat$month))

mdat = mdat %>% arrange(month)

mdat$month_norm = 2 * pi * mdat$month / max(mdat$month,na.rm=T)
mdat$month_cyclical = cos(mdat$month_norm)

data1 = data %>% filter(grepl('SL',sample_id))
merged1 = inner_join(data1,mdat,by=c('sample_id'='ha_sample_id'))%>% distinct

data2 = data %>% filter(!grepl('SL',sample_id))# %>% mutate(sample_id= gsub('\\.','-',sample_id))
merged2 = inner_join(data2,mdat,by=c('sample_id'='sample_id')) %>% distinct

merged1 = bind_rows(merged1,merged2) %>% filter(site_name_description!='n/a')

merged1 = inner_join(merged1,readcountspostqc %>% filter(V2!=0)) 

data_wide = merged1 %>% mutate(val=1) %>% dcast(sample_id + V2 + LOCATION_TYPE ~ GENE,value.var= 'val') 

allids = readRDS('~/Dropbox (Mason Lab)/radx/data_packet_feb23/bulk_mdat_matchedids.rds') %>% filter(!(sample_id %in% data_wide$sample_id))
allids$LOCATION_TYPE[allids$LOCATION_TYPE == 'WASTEWATER TREATMENT PLANT'] = "Wastewater Treatment Plant"
allids$LOCATION_TYPE[allids$LOCATION_TYPE == 'HOSPITAL'] = "Hospital/Medical Campus"
allids$LOCATION_TYPE[allids$LOCATION_TYPE == 'CAMPUS BASIN'] = "University Campus"
allids$LOCATION_TYPE[allids$LOCATION_TYPE == 'DORMITORY'] = "Dormitory"
allids$LOCATION_TYPE[allids$LOCATION_TYPE == 'CAMPUS BASIN / DORMITORIES'] = "University Campus"

data_wide = bind_rows(data_wide,allids)
data_wide[is.na(data_wide)] = 0

# heatmap

# filter out VFs with a prev of 1 

locvec = data_wide %>% arrange(LOCATION_TYPE)  %>% select(LOCATION_TYPE) %>% unlist %>% unname
data_wide_hm = data_wide %>% column_to_rownames('sample_id') %>% arrange(LOCATION_TYPE)%>% select(-LOCATION_TYPE,-SEWERSHED_SCALE) 

tokeep= colSums(data_wide_hm %>% select(-V2)) %>% data.frame %>% filter(.>1) %>% rownames
data_wide_hm = data_wide_hm %>% select(all_of(tokeep))

data_wide_anno = data_wide %>% select(sample_id,LOCATION_TYPE,SEWERSHED_SCALE,V2)
temp = rowSums(data_wide_hm) %>% data.frame %>% rownames_to_column('sample_id') %>% dplyr::rename(`VF Count` = ".")
data_wide_anno = inner_join(data_wide_anno,temp)%>% column_to_rownames('sample_id') %>% mutate(`VFs/10M reads` = `VF Count`/(V2/10000000))
data_wide_anno=data_wide_anno[rownames(data_wide_hm),]
row_ha = rowAnnotation(`VFs/10M reads` = data_wide_anno$`VF Count`,LOCATION_TYPE = data_wide_anno$LOCATION_TYPE)

#res2 = res %>% t %>% data.frame(check.names=F) %>%  select(all_of(colnames(data_wide_hm)))

col_fun  = circlize::colorRamp2(c(0,1,5), c('black','blue','darkred'))
a = Heatmap(data_wide_hm,show_column_names = T,show_row_names = F,left_annotation = row_ha,split = locvec,row_title_rot = 0,col=col_fun,height=8)
colors = structure(1:2, names = c("0","1"))
#b = Heatmap(res2,show_column_names = T,show_row_names = T,col=colors,height=4)

pdf('~/Dropbox (Mason Lab)/radx/plots/VF_heatmap_vfdb.pdf',width=14,height=10)
a #%v% b
dev.off()

plae <- c( `University Campus`= "#E41A1C", `Dormitory` = "#377EB8", "Hospital/Medical Campus" =  "#984EA3", `Wastewater Treatment Plant` = "#FF7F00")

### PREVALENCE BY GROUP
temp = merged1 %>% select(sample_id,LOCATION_TYPE) %>% distinct
mdatsub = mdat %>% filter(sample_id %in% allids$sample_id | ha_sample_id %in% allids$sample_id)  %>% select(sample_id,ha_sample_id,LOCATION_TYPE) %>% distinct
a = mdatsub%>% filter(grepl('SL',ha_sample_id))  %>% select(ha_sample_id,LOCATION_TYPE) %>% dplyr::rename(sample_id = ha_sample_id)
b = mdatsub%>% filter(!grepl('SL',ha_sample_id)) %>% select(-ha_sample_id)

temp = bind_rows(temp,a,b) %>% distinct 
data_wide_prev = inner_join(data_wide %>% melt(id.vars = colnames(data_wide %>% select(V2,sample_id,LOCATION_TYPE,SEWERSHED_SCALE))),temp)
prev = data_wide_prev %>% group_by(LOCATION_TYPE,variable) %>% summarise(Freq = sum(value)) %>% arrange(desc(Freq))%>% filter(Freq>=3)  
prev2 = prev %>% group_by(variable) %>% summarise(Freq=sum(Freq)) %>% arrange(desc(Freq)) %>% mutate(variable=as.character(variable))
prev$variable=factor(prev$variable,levels=(prev2$variable))
prev = inner_join(prev,prev %>% ungroup %>% group_by(variable) %>% mutate(total = sum(Freq)) %>% select(-Freq)) %>% mutate(prop = Freq/total)
plot = ggplot(prev,aes(x = prop,y=variable,fill=LOCATION_TYPE)) + geom_bar(stat='identity',width=.9) + theme_bw() + ylab('') + xlab('Prevalence') + ylab('') + scale_fill_manual(values=plae) + geom_text(aes(label = total, x = 1 + 0.1,size=8, y = variable), size = 3, color = "black") + coord_flip() + theme(legend.position = 'none',axis.text.x = element_text(angle=45,hjust =1))

ggsave(plot=plot,'~/Dropbox (Mason Lab)/radx/plots/VF_prev_vfdb.pdf',width=5,height=2)

# most sample type with genes 
data_wide_sum = data_wide %>% select(-SEWERSHED_SCALE) %>% melt(id.vars = colnames(data_wide %>% select(V2,sample_id,LOCATION_TYPE))) %>% group_by(sample_id,LOCATION_TYPE,V2) %>% summarise(count  = sum(value)) %>% mutate(prop = count/(V2/10000000)) %>% mutate(prop = if_else(is.na(prop),0,prop))
data_wide_sum$LOCATION_TYPE = factor(data_wide_sum$LOCATION_TYPE,levels = c('Hospital/Medical Campus','Dormitory','University Campus','Wastewater Treatment Plant'))
ggplot(data_wide_sum,aes(x = LOCATION_TYPE,y=prop)) + geom_quasirandom(size=3) + theme_bw() + xlab('') + ggtitle('VF identification by site') + ylab('VFs/10M reads') + stat_compare_means(method = 'wilcox.test',comparisons = list(c('Hospital/Medical Campus','Dormitory'),c('Dormitory','University Campus'),c('Dormitory','Wastewater Treatment Plant'))) + theme(axis.text.x = element_text(angle = 45,hjust=1))
ggsave('~/Dropbox (Mason Lab)/radx/plots/locrich_vfdb.pdf',width=3.5,height=8)
dev.off()

# VF over time averaged

temp = merged1 %>% select(sample_id,month2,month) %>% distinct

mdatsub = mdat %>% filter(sample_id %in% allids$sample_id | ha_sample_id %in% allids$sample_id)  %>% select(sample_id,ha_sample_id,month,month2) %>% distinct
a = mdatsub%>% filter(grepl('SL',ha_sample_id))  %>% select(ha_sample_id,month,month2) %>% dplyr::rename(sample_id = ha_sample_id)
b = mdatsub%>% filter(!grepl('SL',ha_sample_id)) %>% select(-ha_sample_id)

temp = bind_rows(temp,a,b) %>% distinct 

data_wide_sum_month = inner_join(data_wide_sum,temp) %>% mutate(prop = count/(V2/10000000))

ggplot(data_wide_sum_month,aes(x = month2,y=count,color=LOCATION_TYPE))+ geom_quasirandom(size=3) + theme_bw() + ylab('VFs/10M reads') + ggtitle('VF richness over time')+xlab('') + scale_color_manual(values=plae) + theme(legend.position ="None")
ggsave('~/Dropbox (Mason Lab)/radx/plots/VFs_over_time_vfdb.pdf',width=5,height=4)


