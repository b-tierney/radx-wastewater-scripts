# predictability in clinic

#### LOAD IN THE CODE USED TO MAKE THE BIG PLOT
#### LINE UP THE DATES ON THIS PROPORTION PLOT AND THE RAW COUNT PLOT, OR MAYBE ADD LINES FOR RAW COUNTS
#### ADD A LINE FOR CASES
#### ANNOTATE WHEN DELTA FIRST IDENTIFIED
#### ANNOTATE WHEN DELTA FIRST IDENTIFIED IN WASTEWATER BY LOCATION

library(tidyverse)

# load the yamina data
setwd('~/Dropbox (Mason Lab)/radx/data_packet_feb23/')

data = read.delim('~/Dropbox (Mason Lab)/radx/data_packet_feb23/yamina/yamina_data4.csv',sep=',',header=T) %>% filter(!is.na(Location)) #%>% mutate(NGS = if_else(NGS == '',qPCR.Screen.Interpretation,NGS))#%>% filter(Location == 'UM-CG')
data$Collection.Date = as.Date(data$Collection.Date,"%m/%d/%y") 
#data$N = as.numeric(data$N)

vars = read.delim('~/Dropbox (Mason Lab)/radx/data_packet_feb23/yamina/variant_mapping.csv',sep=',',header=F) 
  earldates = read.csv('~/Dropbox (Mason Lab)/radx/data_packet_feb23/yamina/earliest_dates_variants_observed.csv') %>% mutate(date_observed = as.Date(date_observed,"%m/%d/%y"))
colnames(vars) = c('Lineage','var')
data = left_join(data,vars %>% distinct,by='Lineage')
#data = data %>% mutate(var = if_else(NGS == 'SGTA' | NGS== "SGTF",'Ambiguous',NGS))

data = data %>% dplyr::rename(date = Collection.Date,location = Location, lineage = var) %>% select(date,location,lineage,Sampletype)

# load the foox data
lineagetrack = read.delim('~/Dropbox (Mason Lab)/radx/data_packet_feb23/ARTIC/freyja_outs_full.csv',header=T,sep=',')
df_meta = read.delim("ARTIC/foox_apr2023_all_plates_metadata.bb.tsv", sep="\t")

df_meta$date = as.Date(df_meta$date,"%Y-%m-%d")

lineagetrack = inner_join(lineagetrack,df_meta) %>% select(date,location,lineage,fraction) %>% rename(n = fraction) %>% group_by(date,lineage) %>% summarise(n = n()) %>% mutate(n = n / sum(n)) %>% mutate(Sampletype = 'Wastewater Tracking')

data2 = data %>% group_by(date,lineage,Sampletype) %>% summarise(n = n()) %>% mutate(n = n / sum(n))

toplot = bind_rows(lineagetrack,data2)#%>% filter(date > "2021-01-01")
toplot$lineage=factor(toplot$lineage,levels = c('A','Alpha','Beta','Epsilon','Iota','Gamma','Zeta','Eta','Delta','Mu','BA.1','BA.2','BA.4or5','Ambiguous'))

lineage.colors = c("Alpha"="#3F4CCB",
                   "Beta"="#492AB5",
                   "Gamma"="#4271CE",
                   "Delta"="#4C8FC0",
                   "BA.1"="#DB2823",
                   "BA.1or2"="#bd1915",
                   "BA.2"="#9e110d",
                   "BA.4or5" = "#630300",
                   "BA.5"= "#1c0100",
                   "Lambda"="#E68033",
                   "Mu"="#ff5500",
                   "Epsilon"="#A0BE59",
                   "Eta"="#BBBC49",
                   "Iota"="#E19F3A",
                   "Kappa"="#85BA6F",
                   "Ambiguous"="#999999")


toplot2 = left_join(toplot,earldates) %>% mutate(lineage = if_else(is.na(lineage),'Ambiguous',lineage))

toplot2 = toplot2 %>% filter(is.na(date_observed) | !is.na(date_observed) & date_observed<date)


ints = toplot2 %>% ungroup%>% select(lineage,date_observed) %>% distinct %>% arrange()


ggplot(toplot2,aes(x = date,y= n,fill=lineage)) + geom_area(size=5,alpha=1) + theme_cowplot() + facet_grid(Sampletype ~ .,scales='free_y',space = 'free_x') + ylab('Cumulative VoC proportion') + xlab('') + scale_fill_manual(values=lineage.colors) + theme(legend.position = 'none') 
ggsave('~/Dropbox (Mason Lab)/radx/plots/clinical_variant_map.pdf',width=5,height=5)


#data2 = full_join(data2 %>% group_by(lineage,Sampletype) %>% slice_min(date,n=1) %>% select(lineage,date,Sampletype) %>% filter(Sampletype == 'Patient')%>% dplyr::rename(Patient = date) %>% ungroup %>% select(-Sampletype),data2 %>% group_by(lineage,Sampletype) %>% slice_min(date,n=1) %>% select(lineage,date,Sampletype) %>% filter(Sampletype == 'Student') %>% dplyr::rename(Student = date) %>% ungroup %>% select(-Sampletype))


#### get date first detected by location
#diffplot = inner_join(inner_join(lineagetrack,earldates) %>% group_by(lineage) %>% filter(date_observed<date)  %>% distinct %>% group_by(lineage) %>% slice_min(date,n=1) %>% slice_max(n,n=1) %>% select(lineage,date) %>% rename(Wastewater = date),data2) %>% filter(lineage != 'Other') %>% mutate(diff_p = Wastewater - Patient,diff_s = Wastewater - Student) %>% arrange(Patient) %>% filter(lineage %in% c('Delta','BA.1','BA.2'))

#ggplot(diffplot,aes(x = diff, y= reorder(lineage,Patient))) + geom_point(size=3) + xlim(-60,60) + theme_cowplot() + geom_vline(linetype = 'dashed',color='red',xintercept = 0) + geom_hline(aes(yintercept = reorder(lineage,Patient)),linetype = 'dashed') + ylab('') + ggtitle('Day variant detected in wastewater relative to clinic') + xlab('Day offset')
#ggsave('~/Dropbox (Mason Lab)/radx/plots/clinical_vs_wawa_timing.pdf',width=7,height=5)








