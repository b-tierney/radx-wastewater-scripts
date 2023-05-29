---
title: 'UMiami Wastewater ARTIC Total Data Viz'
author: "Jonathan Foox"
date: "April 18, 2023"
output:
  html_document:
    df_print: paged
  rmdformats::html_docco:
    gallery: yes
    lightbox: yes
    thumbnails: no
---

```{=html}
<style type="text/css">
.main-container {
  max-width: 1800px;
  margin-left: auto;
  margin-right: auto;
}
</style>
```

```{r, message=F}
library(tidyr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(cowplot)
library(ggpubr)
options(scipen=999)
library(scales)
library(ComplexHeatmap)
library(circlize)
library(viridis)
```

```{r, message=F}
source("/home/cem2009/R/cem_helper.R")
source("/home/cem2009/R/cem_plot.R")
source("/home/cem2009/R/cem_plot2.R")

theme_jon = theme_linedraw() + 
  theme(strip.background=element_rect(fill="grey80", colour="grey50", size=0.2), strip.text.x=element_text(colour="black"), strip.text.y=element_text(colour="black")) +
  theme(legend.position="bottom", 
        legend.title=element_text(size=11), 
        legend.text=element_text(size=10), 
        axis.text.x = element_text(size=10, angle=90, hjust=1, vjust=0.5), 
        axis.text.y=element_text(size=10), 
        axis.title=element_text(size=12), 
        axis.title.y=element_text(vjust=1), 
        plot.title = element_text(size=12, hjust = 0.5), 
        strip.text = element_text(size=12), 
        panel.grid.major = element_line(colour = "grey98"), 
        panel.grid.minor = element_blank())

# see freyja section
#lineage.colors = c("Alpha"="#3F4CCB", "Beta"="#492AB5", "Gamma"="#4271CE", "Delta"="#4C8FC0", "Omicron"="#DB2823", "Lambda"="#E68033", "Mu"="#E2562B", "Epsilon"="#A0BE59", "Eta"="#BBBC49", "Iota"="#E19F3A", "Kappa"="#85BA6F", "Ambiguous"="#999999")#, "A"="#8b99b3", "Theta"="#88b891", "Zeta"="#b38ba0")
```

# Metadata table to convert sample names

```{r}
df_meta = read.csv("foox_apr2023_all_plates_metadata.bb.tsv", sep="\t")

# deal with dups in a silly way for now
df_meta$rep = "1"
df_meta$sample2 = paste0(df_meta$date, "__", df_meta$location, "__", df_meta$rep)
df_meta[duplicated(df_meta$sample2), ]$rep = "2"
df_meta$sample2 = paste0(df_meta$date, "__", df_meta$location, "__", df_meta$rep)
df_meta[duplicated(df_meta$sample2), ]$rep = "3"
df_meta$sample2 = paste0(df_meta$date, "__", df_meta$location, "__", df_meta$rep)

# add year/month col
df_meta$date = as.Date(df_meta$date)
df_meta$date2 = df_meta$date
df_meta$date2 = gsub("-\\d\\d$", "", df_meta$date2)
df_meta
```

```{r}
df_meta_sites = read.csv("site_table.tsv", sep="\t")
```

# Kraken Taxonomic Breakdown

```{r, message=F}
df_kraken = read.csv("tables/kraken__taxonomy.csv")
df_kraken$other_viruses = df_kraken$viruses - df_kraken$sars_cov2

taxaCol = c("sars_cov2" = "#ae017e", "human" = "#6797CB", "fungi" = "#53BA98", "bacteria" = "#B77350", "archaea" = "#B09C42",  "other_viruses" = "#C5688E", "not in database"="#999999")
df_kraken2 = df_kraken[, c("sample", "archaea", "bacteria", "fungi", "human", "sars_cov2", "other_viruses", "unclassified")]
df_kraken2 = setRemoveRownames(df_kraken2, 1)
df_kraken2scaled = df_kraken2
df_kraken2scaled = df_kraken2scaled / rep.col( rowSums(df_kraken2scaled), ncol(df_kraken2scaled)-1)

myOrder = getTSPOrder(df_kraken2scaled)

df_kraken3 = addRownameColumn(df_kraken2[myOrder, ], "sample")
df_kraken4 = melt(df_kraken3, id.vars=c("sample"), measure.vars=c("bacteria","archaea","fungi","human","sars_cov2","other_viruses","unclassified"), variable.name="Classification", value.name="Value")
df_kraken4$sample = factor(df_kraken4$sample, levels=myOrder)
df_kraken4$sampleType = ifelse(grepl("-VarSkip", df_kraken4$sample), "VarSkip", "V3")

# change unclassified
df_kraken4$Classification = gsub("unclassified", "not in database", df_kraken4$Classification)
```

```{r, fig.width=11, fig.height=6}

taxa_plot_num_upper = ggplot(df_kraken4, aes(x=sample, y=Value, fill=Classification)) + 
  geom_bar(stat="identity", position="stack") + 
  geom_hline(yintercept = 0) +
  theme_jon + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(1, 1, 0, 0), "cm"),
        legend.position="top") +
  scale_fill_manual(values=taxaCol) + scale_color_manual(values=taxaCol) +
  scale_y_continuous(labels=comma) +
  coord_trans(ylim=c(5000000,20000000)) +
  xlab("") + ylab("")

taxa_plot_num_lower = ggplot(df_kraken4, aes(x=sample, y=Value, fill=Classification)) + 
  geom_bar(stat="identity", position="stack", show.legend=F) + 
  geom_hline(yintercept = 0) +
  theme_jon + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0, 1, 0, 1), "cm")) +
  scale_fill_manual(values=taxaCol) + scale_color_manual(values=taxaCol) +
  scale_y_continuous(labels=comma) +
  coord_trans(ylim=c(0,5000000)) +
  xlab("") + ylab("# Reads\n")

taxa_plot_fct = ggplot(df_kraken4, aes(x=sample, y=Value, fill=Classification)) + 
  geom_bar(stat="identity", position="fill") + 
  geom_hline(yintercept = 0) +
  theme_jon + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0, 1, 0, 1), "cm"),
        legend.position="none") +
  scale_fill_manual(values=taxaCol) + scale_color_manual(values=taxaCol) +
  xlab("") + ylab("Fraction\n") + 
  scale_y_continuous(breaks=c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1)) + coord_trans(y="sqrt")

plot_grid(taxa_plot_num_upper, taxa_plot_num_lower, taxa_plot_fct, ncol=1, align = "v", rel_heights = c(1, 1.1, 1))
ggsave(file="figures/2023-04-17__kraken__sorted.png", width=11, height=8)
```

```{r}
df_kraken_dateOrder = df_kraken
df_kraken_dateOrder = merge(df_kraken_dateOrder, df_meta[,c("sample","sample2")])
df_kraken_dateOrder = df_kraken_dateOrder[,-1]
df_kraken_dateOrder = df_kraken_dateOrder %>% select(sample2, everything())
df_kraken_dateOrder = df_kraken_dateOrder[, c("sample2", "archaea", "bacteria", "fungi", "human", "sars_cov2", "other_viruses", "unclassified")]


df_kraken_dateOrder2 = melt(df_kraken_dateOrder, id.vars=c("sample2"), 
                            measure.vars=c("bacteria","archaea","fungi","human","sars_cov2","other_viruses","unclassified"), 
                            variable.name="Classification", value.name="Value")

df_kraken_dateOrder2$sample2 = factor(df_kraken_dateOrder2$sample2, levels=unique(sort(df_kraken_dateOrder2$sample2)))

# change unclassified
df_kraken_dateOrder2$Classification = gsub("unclassified", "not in database", df_kraken_dateOrder2$Classification)

```

```{r, fig.width=11, fig.height=6}

taxa_plot_num_upper = ggplot(df_kraken_dateOrder2, aes(x=sample2, y=Value, fill=Classification)) + 
  geom_bar(stat="identity", position="stack") + 
  geom_hline(yintercept = 0) +
  theme_jon + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(1, 1, 0, 0), "cm"),
        legend.position="top") +
  scale_fill_manual(values=taxaCol) + scale_color_manual(values=taxaCol) +
  scale_y_continuous(labels=comma) +
  coord_trans(ylim=c(5000000,20000000)) +
  xlab("") + ylab("")

taxa_plot_num_lower = ggplot(df_kraken_dateOrder2, aes(x=sample2, y=Value, fill=Classification)) + 
  geom_bar(stat="identity", position="stack", show.legend=F) + 
  geom_hline(yintercept = 0) +
  theme_jon + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0, 1, 0, 1), "cm")) +
  scale_fill_manual(values=taxaCol) + scale_color_manual(values=taxaCol) +
  scale_y_continuous(labels=comma) +
  coord_trans(ylim=c(0,5000000)) +
  xlab("") + ylab("# Reads\n")

taxa_plot_fct = ggplot(df_kraken_dateOrder2, aes(x=sample2, y=Value, fill=Classification)) + 
  geom_bar(stat="identity", position="fill") + 
  geom_hline(yintercept = 0) +
  theme_jon + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0, 1, 0, 1), "cm"),
        legend.position="none") +
  scale_fill_manual(values=taxaCol) + scale_color_manual(values=taxaCol) +
  xlab("") + ylab("Fraction\n") + 
  scale_y_continuous(breaks=c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1)) + coord_trans(y="sqrt")

plot_grid(taxa_plot_num_upper, taxa_plot_num_lower, taxa_plot_fct, ncol=1, align = "v", rel_heights = c(1, 1.1, 1))
ggsave(file="figures/2023-04-17__kraken__timeorder.png", width=11, height=8)
```

# Genome Coverage (mosdepth)



```{r import coverage table}
df_cov = read.csv("tables/mosdepth_allsamples.tsv", sep="\t")
df_cov = dcast(df_cov, sample~amplicon, value.var="coverage")

# get rows sorted by date
df_cov = merge(df_cov, df_meta[,c("sample","sample2")], by="sample")
df_cov = df_cov[, -c(1)]
df_cov = df_cov %>% select(sample2, everything())
df_cov = df_cov[order(df_cov$sample2), ]

rownames(df_cov) = df_cov$sample2


# reduce just to WW samples
#df_cov = subset(df_cov, grepl("_W", df_cov$sample))

# sort lexicographically
df_cov = df_cov[, stringr::str_sort(colnames(df_cov[,-1]), numeric = T)]

# turn into matrix
df_covm = as.matrix(df_cov)
df_covm[1:5, 1:5]
```

```{r}
# add new object to split genomecov map into good/bad coverage
df_cov_goodbad = data.frame(read.csv("samples_genomecov_min75pct.txt", sep="\t", header=F))
colnames(df_cov_goodbad) = "sample"
passSamples = df_cov_goodbad$sample

df_cov_goodbad = data.frame(dcast(read.csv("tables/mosdepth_allsamples.tsv", sep="\t"), sample~amplicon, value.var="coverage")$sample)
colnames(df_cov_goodbad) = "sample"
df_cov_goodbad$passfail = ifelse(df_cov_goodbad$sample %in% passSamples, "pass", "fail")
df_cov_goodbad$passfail = factor(df_cov_goodbad$passfail, levels=c("pass", "fail"))
df_cov_goodbad
```

```{r}
df_cov_pct_amplicons = read.csv("tables/amplicon_min10xcov.csv", header=F)
colnames(df_cov_pct_amplicons) = c("sample", "num_amplicons")
df_cov_pct_amplicons = merge(df_cov_pct_amplicons, df_meta[,c("sample","sample2")])
df_cov_pct_amplicons$num_amplicons_pct = df_cov_pct_amplicons$num_amplicons / 98
df_cov_pct_amplicons
```

```{r import metadata on genes and batches}
# metadata

### gene colors/labels
gene.colors = c("Spike"="#de2d26", "Nuc"="#c51b8a", "Env"="#31a354", "Mem"="#3182bd", "ORF1ab"="#ffffe5", "ORF3a"="#fff7bc", "ORF67ab"="#fee391", "ORF8"="#fe9929", "ORF10"="#ec7014")
gene_split = read.csv("../artic-nCov-2019.v3.insert.proper.renamedchr.bed", sep="\t", header=F)
gene_split = gene_split[, c(4,7)]
colnames(gene_split) = c("amplicon", "gene")
gene_split$gene = factor(gene_split$gene, levels=c( "ORF1ab", "Spike", "ORF3a", "Env", "Mem", "ORF67ab", "ORF8", "Nuc", "ORF10"))
gene_split
```


```{r plot heatmap, fig.width=10, fig.height=8}

hm_cov = Heatmap(log10(df_covm+1),
        name="Cov", show_heatmap_legend = T,
        col = viridis(50),
        cluster_rows = F,
        cluster_row_slices = F,
        show_row_names = F,
        #row_split = df_cov_goodbad$passfail,
        
        row_names_gp = gpar(fontsize=7),
        row_title_gp = gpar(fontsize=10),
        row_title_rot = 0,
        
        cluster_columns = F,
        show_column_names = T,
        column_names_gp = gpar(fontsize=6),
        column_title_gp = gpar(fontsize=8),
        column_title_rot = 90,
        column_split = gene_split$gene,
        #rect_gp = gpar(col = "white", lwd = 0.1),
        column_gap = unit(2, "mm"), border = TRUE,
        right_annotation = rowAnnotation(pct  = anno_barplot(df_cov_pct_amplicons$num_amplicons_pct, gp = gpar(col = ifelse(df_cov_pct_amplicons$num_amplicons_pct > 0.75, "#3288bd", "#d53e4f")))))
        


draw(hm_cov, padding = unit(c(10, 1, 1, 1), "mm"))

pdf(file="figures/2023-04-17__genomecov__timeorder.pdf", width=10, height=6)
draw(hm_cov, padding = unit(c(10, 1, 1, 1), "mm"))
decorate_annotation("pct", {   grid.lines(x = unit(c(0.75, 0.75), "npc"), y = unit(c(0, 855), "native"), gp = gpar(lty = 2, lwd = 1.5, col = "yellow")) })
dev.off()
```

# Variant Viz

```{r VOC metadata}
vocs = read.csv("../FORPUB/COVARIANTS_sig_mutations.withGeneName.withBA2-5.csv")
vocs$gene = factor(vocs$gene, levels = c("ORF1a", "ORF1b", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF8", "N", "ORF9b"))

# separate rows so that uniquetons are shown first, then multiples
vocs$dupes = duplicated(vocs$AAsub) | duplicated(vocs$AAsub, fromLast=TRUE)

#order_of_VOCs = c("omicron", "delta", "alpha", "beta", "gamma", "kappa", "eta", "iota", "lambda", "mu", "epsilon")
#vocs$name = factor(vocs$name, levels = order_of_VOCs)
vocs = vocs[order(vocs$dupes, vocs$name), ]
vocs
```

```{r VOC annotations}
vocs2 = dcast(vocs, AAsub~name, length)
vocs2
```

```{r get variant info for all samples}
df_muts_all = data.frame()

for ( wwsample0 in list.files(path="/athena/masonlab/scratch/projects/metagenomics/covid/analysis/miami_wastewater/latest_all_plates/varcall/outs", 
                              include.dirs=F, recursive = F, pattern=glob2rx("ivarAndLoFreq_VEP_*"))) {
  wwsample = gsub("ivarAndLoFreq_VEP_", "", wwsample0)
  wwsample = gsub(".csv", "", wwsample)
  
  vepfile = paste0("varcall/outs/",wwsample0)
  df_vep  = read.table(vepfile, sep=",", header=T)

  ### replace "-"s with NAs
  df_vep = na_if(df_vep, "-")
  
  
  
  ### transform into matrix with only the data we need
  df_sample_muts = data.frame(gene_mut_loc = str_split_fixed(df_vep$Location,"[:-]+", n=3),
                         prot_mut_loc = df_vep$Protein_position, 
                         AAs = str_split_fixed(df_vep$Amino_acids, "/",2), 
                         Conseq = df_vep$Consequence,
                         genes = df_vep$SYMBOL)
  df_sample_muts = na_if(df_sample_muts,'')
  
  ### remove rows with no protein position value
  df_sample_muts = df_sample_muts %>% filter(!grepl("^-", prot_mut_loc))
  
  ### adjust AA labels and remove non-protein coding sites
  df_sample_muts$AAs.2 = ifelse(df_sample_muts$Conseq == "synonymous_variant", df_sample_muts$AAs.1, df_sample_muts$AAs.2)
  df_sample_muts$AA_mut = paste(df_sample_muts$AAs.1,df_sample_muts$prot_mut_loc,df_sample_muts$AAs.2, sep = "")
  df_sample_muts = subset(df_sample_muts, AA_mut != "NANANA")
  
  ### adjust further for deletions, which are represented as a block, but should be broken up to one AA at a time
  write.csv(df_sample_muts, file="df_sample_muts.tmp", row.names = F, quote = F)
  system("python ../run_split_deletion_mutations_for_Rmd.py")
  df_sample_muts = read.csv("df_sample_muts.updated.tmp")
  
  
  
  if(nrow(df_sample_muts) > 0 ) {
    ### append VOC label(s) to table
    df_sample_muts = dplyr::left_join(df_sample_muts, vocs[, c("name", "AAsub")], by = c('AA_mut'='AAsub'))
    
    ### get FREQUENCY information from parsed IVAR outputs
    ivarfile = paste0("varcall/outs/ivarAndLoFreq_VCF_",wwsample,".csv")
    df_ivar = read.csv(ivarfile, header=F)
    colnames(df_ivar) = c("POS", "REF", "ALT", "TOTAL_DP", "ALT_FREQ")
    
    # account for 0-based and 1-based differences
    df_ivar$POS = ifelse(grepl("-", df_ivar$ALT), df_ivar$POS+1, df_ivar$POS)
    
    df_ivar$gene_mut = paste0(df_ivar$REF, df_ivar$POS, df_ivar$ALT)
    df_ivar = df_ivar[, c("POS", "gene_mut", "ALT_FREQ", "TOTAL_DP")]
    colnames(df_ivar) = c("gene_pos", "gene_mut", "freq", "cov")
    df_ivar$gene_pos = as.character(df_ivar$gene_pos)
      
    
    ### combine with sample mutation table
    df_sample_muts$gene_mut_loc.2 = as.character(df_sample_muts$gene_mut_loc.2)
    df_sample_muts = left_join(df_ivar, df_sample_muts, by=c("gene_pos"="gene_mut_loc.2"), copy=T)
    df_sample_muts = unique(df_sample_muts)
    
    ### filter sample mutations into VOC matches
    df_matches = df_sample_muts %>% filter(!is.na(name))
    
    if(nrow(df_matches) > 0 ) {
      ### add sample name
      df_matches$sample = wwsample
      df_muts_all = rbind(df_matches, df_muts_all)    
    }  
  }
  
}
df_muts_all
```

```{r transform into matrix for heatmap}
# convert tall to wide
df_muts_all2 = dcast(unique(df_muts_all[,-c(2,13)]), sample~AA_mut, value.var="freq", max)


# replace -Infs
df_muts_all2[mapply(is.infinite, df_muts_all2)] <- NA

# fill in NAs based on whether that position was covered or not
# 1. read in genomecov as object
# 2. for each sample, for each AA column, look for an NA
# 3. consult AA->nuc table (derived from unique df_muts_all nuc/AA object)
# 4. look up that nuc in the genomecov object
# 5. if >0, then 0.0, elif=0 then -1 or whatever
AAtoNuc = unique(df_muts_all[, c("gene_pos", "AA_mut")])

for(samplename in df_muts_all2$sample) {
  currentgenomecov = read.csv(paste0("bwa/outs/",samplename,".ivar.genomecov"), sep="\t", header=F)
  colnames(currentgenomecov) = c("genome", "base", "cov")
  for(AAmutname in colnames(df_muts_all2)[2:ncol(df_muts_all2)]) {
    if(is.na(df_muts_all2[df_muts_all2$sample==samplename, AAmutname])) {
      nucpos = as.integer(AAtoNuc[AAtoNuc$AA_mut==AAmutname, ]$gene_pos[1]) # in case it shows up twice
      coverage = as.integer(currentgenomecov[currentgenomecov$base==nucpos, ]$cov)
      df_muts_all2[df_muts_all2$sample==samplename, AAmutname] = ifelse(coverage > 9, 0, -.01)
    }
  }
}
```

```{r}
df_muts_all2m = df_muts_all2
df_muts_all2m = merge(df_muts_all2m, df_meta[, c("sample", "sample2")])
df_muts_all2m = df_muts_all2m %>% select(sample, sample2, everything())

#df_muts_all2m = separate(df_muts_all2m, col="sample", sep="-", into=c("site", "date", "extra"), remove=T)
#df_muts_all2m$sample2 = df_muts_all2m$sample#paste0(df_muts_all2m$date, "-", df_muts_all2m$site, "-", df_muts_all2m$extra)

df_muts_all2m = df_muts_all2m[order(df_muts_all2m$sample2), ]


#df_muts_all2m = df_muts_all2m[,-c(1,2,3)]
#df_muts_all2m = df_muts_all2m %>% select(sample, everything())

#back to normal
rownames(df_muts_all2m) = df_muts_all2m[,c("sample2")]
df_muts_all2m = as.matrix(df_muts_all2m[,-c(1,2)])
df_muts_all2m[1:5, 1:5]
rownames(df_muts_all2m)
```

```{r subset VOCs to variants observed}
vocs3 = subset(vocs2, AAsub %in% colnames(df_muts_all2))
vocs3$duplicated = ifelse(rowSums(vocs3[,2:ncol(vocs3)])>1, "Recurrent", "Unique")
vocs3$duplicated = factor(vocs3$duplicated, levels=c("Unique", "Recurrent"))
vocs3


# account for different deltas
delta_21A = c("T19R", "E156-", "F157-", "R158G", "L452R", "T478K", "D614G", "P681R", "D950N", "P314L", "G662S", "P1000L", "I82T", "D63G", "R203M", "D377Y", "S26L", "V82A", "T120I", "D119-", "F120-", "T60A")
delta_21I = c("A222V", "P1640L", "A3209V", "V3718A", "T3750I")
delta_21J = c("G215C", "A1306S", "V2930L", "T3255I", "T3646A", "A1918V", "T40I")
vocs3$delta = ifelse(vocs3$AAsub %in% delta_21A, 1, ifelse(vocs3$AAsub %in% delta_21I, 2, ifelse(vocs3$AAsub %in% delta_21J, 3, 0)))
```

```{r}
hm_left_annotation = rowAnnotation(Alpha        = vocs3$alpha, 
                                   Beta         = vocs3$beta, 
                                   Delta        = vocs3$delta,
                                   Epsilon      = vocs3$epsilon,
                                   Eta          = vocs3$eta, 
                                   Gamma        = vocs3$gamma,
                                   Iota         = vocs3$iota, 
                                   Kappa        = vocs3$kappa, 
                                   Lambda       = vocs3$lambda,
                                   Mu           = vocs3$mu,
                                   OmicronBA1   = vocs3$omicronBA1,
                                   OmicronBA2   = vocs3$omicronBA2,
                                   OmicronBA212 = vocs3$`omicronBA2.12.1`,
                                   OmicronBA4   = vocs3$omicronBA4,
                                   OmicronBA5   = vocs3$omicronBA5,
                                   
                                   col = list(Alpha = c("0"="#ffffff", "1"="#5d85e6"),
                                              Gamma = c("0"="#ffffff", "1"="#2d34bd"),
                                              Lambda = c("0"="#ffffff", "1"="#f6e8c3"),
                                              Beta = c("0"="#ffffff", "1"="#d16ef5"),
                                              Epsilon = c("0"="#ffffff", "1"="#ebce3b"),
                                              Iota = c("0"="#ffffff", "1"="#c9af28"),
                                              Eta = c("0"="#ffffff", "1"="#b2df8a"),
                                              Kappa = c("0"="#ffffff", "1"="#33a02c"),
                                              Delta = c("0"="#ffffff", "1"="#da372c", "2"="#ba190f", "3"="#800800"),
                                              Mu = c("0"="#ffffff", "1"="#f2a446"),
                                              OmicronBA1 = c("0"="#ffffff", "1"="#520007"),
                                              OmicronBA2 = c("0"="#ffffff", "1"="#9e110d"),
                                              OmicronBA212 = c("0"="#ffffff", "1"="#9e110d"),
                                              OmicronBA4    = c("0"="#ffffff", "1"="#140002"),
                                              OmicronBA5    = c("0"="#ffffff", "1"="#140002")),
                                           show_legend = F, 
                                   border=T, gap=0, gp = gpar(fontsize=3), 
                                   annotation_name_side = "top", 
                                   annotation_name_gp = gpar(fontsize=10))
```




```{r, fig.width=18, fig.height=12}
hm_var_t_PA = Heatmap(t(df_muts_all2m),
        name="VAF", show_heatmap_legend = T,
        col = colorRamp2(c(-0.01, 0, seq(0.1, 1, 0.1)), c("#ffffff", "#d4ebff", "#FDB863", "#F7A959", "#F29A50", "#ED8B47", "#E77C3D", "#E26D34", "#DD5E2B", "#D74F21", "#D24018", "#CD310F")),
        
        column_split = subset(df_meta, sample2 %in% colnames(t(df_muts_all2m)))$date2,
        #cluster_columns = F,
        column_order = rownames(df_muts_all2m),
        
        cluster_column_slices = F,
        show_column_names = F,
        column_names_gp = gpar(fontsize=10),
        column_title_gp = gpar(fontsize=10),
        column_title_rot = 90,
        show_column_dend = F,
        #column_title = "All Non-Synonymous Variants\n",
        column_gap = unit(1, "mm"), 
        
        

        row_order = unique(subset(vocs, AAsub %in% colnames(df_muts_all2m))$AAsub),
        row_split = vocs3$duplicated,
        row_gap = unit(5, "mm"),
        show_row_names = T,
        row_names_gp = gpar(fontsize=4),
        row_title_gp = gpar(fontsize=10),
        show_row_dend = F,
        row_names_side = "left",
        row_title_rot = 90,
        
        rect_gp = gpar(lwd = 0.2), 
        border = T, 
        left_annotation = hm_left_annotation,
        heatmap_legend_param = list(title = "VAF", at = c(-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1),
        labels = c("Missing Data", "Wild Type", "Var (20%)", "Var (40%)", "Var (60%)", "Var (80%)", "Var (100%)"), border = "black")
        )
draw(hm_var_t_PA, padding = unit(c(10, 10, 10, 10), "mm"))
```

```{r}
pdf(file="figures/2023-04-17__SNPheatmap.pdf", width=18, height=12)
draw(hm_var_t_PA, padding = unit(c(10, 10, 10, 10), "mm"))
dev.off()
```

# Freyja

```{r}
df_freyja = read.csv("tables/freyja_outs.tsv", sep="\t", header=F)
colnames(df_freyja) = c("sample", "VOC", "pct")
df_freyja = df_freyja[complete.cases(df_freyja), ]

# adjust VOC labels
df_freyja$VOC = gsub("Other", "Ambiguous", df_freyja$VOC)

df_freyja = merge(df_freyja, df_meta[,c("sample","sample2")])
df_freyja = separate(df_freyja, col="sample2", sep="__", into=c("date", "site", "rep"), remove=F)
df_freyja = df_freyja[order(df_freyja$sample2), ]

df_freyja$date = as.Date(df_freyja$date)
df_freyja$pct  = as.numeric(df_freyja$pct)

df_freyja = merge(df_freyja, df_meta_sites, by="site")

df_freyja
```

```{r}
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

# #E2562B for mu
```

```{r}
df_freyja = subset(df_freyja, sample != "Babler-13461-037-ARTIC" & sample != "	Babler-13461-038-ARTIC")
```



```{r, fig.width=10, fig.height=6}
gg_stackbar = suppressWarnings(ggplot(subset(df_freyja, site2 != "Marine Campus" & site2 != "Residential Campus")) +
  geom_bar(stat="identity", position="fill", width=5, aes(x=date, y=pct, fill=VOC, group=VOC)) +
  scale_fill_manual(values=lineage.colors) +
  scale_x_date(date_labels="%b %Y", date_breaks="1 month")+
  theme_jon + theme(legend.position = "right", axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  xlab("") + ylab("Fraction") + ggtitle("Wastewater Variant Abundance") +
  theme(strip.text.y=element_text(angle=0)) +
  facet_grid(site2~., space = "free_x", scales = "free_x"))
gg_stackbar
ggsave(file="figures/2023-04-17__freyja__strains-by-aggreg-site.pdf", width=10, height=6)
```

```{r}
df_freyja2 = df_freyja
for(VOC in unique(df_freyja2$VOC)) {
  for(sample in df_freyja2$sample) {
    if( ! any(df_freyja2$sample==sample & df_freyja2$VOC==VOC) ) {
      cursamplesite = df_freyja2[df_freyja2$sample==sample, "site"][1]
      cursample2    = df_freyja2[df_freyja2$sample==sample, "sample2"][1]
      curdate       = df_freyja2[df_freyja2$sample==sample, "date"][1]
      currep        = df_freyja2[df_freyja2$sample==sample, "rep"][1]
      cursite2      = df_freyja2[df_freyja2$sample==sample, "site2"][1]
      df_freyja2[nrow(df_freyja2) + 1,] = list(cursamplesite, sample, VOC, 0, cursample2, curdate, currep, cursite2)
    }
  }
}
```


```{r}
df_freyja2b = subset(df_freyja2, rep=="1")
df_freyja2b = subset(df_freyja2b, site!="WG02" & site!="WG0A" & site!="WG0E" & site !="WG0H" & site !="WG0S" & site !="WG0U")
```


```{r, fig.width=8, fig.height=8}
ggplot(df_freyja2b) +
  geom_area(aes(x=date, y=pct, fill=VOC)) +
  scale_fill_manual(values=lineage.colors) +
  scale_x_date(date_labels="%b %Y", date_breaks="1 month")+
  theme_jon + theme(legend.position = "right", axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  xlab("") + ylab("Fraction") + ggtitle("Wastewater Variant Abundance") +
  theme(strip.text.y=element_text(angle=0), axis.text.y = element_text(size=6)) +
  facet_grid(site~., space = "free_x", scales = "free_x")
ggsave(file="figures/2023-04-17__freyja__strains-by-indiv-sites.pdf", width=8, height=8)
```

```{r}
df_freyja$date2 = format(as.Date(df_freyja$date), "%Y-%m")

```

```{r}
df_freyja2 = reshape2::dcast(df_freyja, date2~VOC, value.var="pct", mean)
df_freyja2[is.na(df_freyja2)] = 0
```

```{r}
df_freyja3 = reshape2::dcast(df_freyja, date~VOC, value.var="pct", mean)
df_freyja3[is.na(df_freyja3)] = 0
```

```{r}
write.csv(df_freyja2, file="tables/permonth_lineage_abundance_table.csv", quote = F, row.names = F)
```

```{r}
df_freycov = read.csv("tables/freyja_outs_with_coverage.tsv", sep="\t", header=F)
colnames(df_freycov) = c("sample", "VOC", "pct")
df_freycov = df_freycov[complete.cases(df_freycov), ]

# adjust VOC labels
df_freycov$VOC = gsub("Other", "Ambiguous", df_freycov$VOC)

df_freycov = merge(df_freycov, df_meta[,c("sample","sample2")])
df_freycov = separate(df_freycov, col="sample2", sep="__", into=c("date", "site", "rep"), remove=F)
df_freycov = df_freycov[order(df_freycov$sample2), ]

df_freycov$date = as.Date(df_freycov$date)
df_freycov$pct  = as.numeric(df_freycov$pct)

df_freycov = merge(df_freycov, df_meta_sites, by="site")

df_freycov
```

```{r}
df_freycov2 = df_freycov
df_freycov2$date2 = format(as.Date(df_freycov2$date), "%Y-%m")
df_freycov2 = reshape2::dcast(df_freycov2, sample+date2~VOC, value.var = "pct")
df_freycov2[is.na(df_freycov2)] = 0


# multiply all by coverage
for(asdf in colnames(df_freycov2)[3:15]){
  df_freycov2[[asdf]] = df_freycov2[[asdf]] * df_freycov2[["coverage"]]
}

# remove coverage
df_freycov2 = df_freycov2[,!grepl("coverage", colnames(df_freycov2))]

# sum up coverages over moth
df_freycov2 = df_freycov2[,c(2:14)] %>% group_by(date2) %>% summarize_all(sum)

# get each row as relative pct
df_freycov2[,c(2:13)] = df_freycov2[,c(2:13)] / rowSums(df_freycov2[,c(2:13)])

df_freycov2
```

```{r}
df_freycov2 %>% group_by(date2) %>% summarise_all(mean)
```

```{r}
rowSums(df_freycov2[,c(3:13)])
```

```{r}
write.csv(df_freycov2, file="tables/permonth_lineage_abundance_table.updated2023-01-10.csv", quote = F, row.names = F)
```
