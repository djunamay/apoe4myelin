##### Plot ACAT1 and 2 in the post-mortem oligodendrocytes

library(SingleCellExperiment)
library("readxl")
library(ggplot2)
library(ggpubr)

# load the data
print('loading the data..')
e = readRDS('../../submission_code_09012021/data/Summary.DE.celltype.rds')
expressed = metadata(e)$expressed.genes
av_expression = readRDS('../../submission_code_09012021/data/Averages.by.celltype.by.individual.rds')

summary = as.data.frame(read_excel('../../submission_code_09012021/data/Overview_PFC_projids_data_384.xlsx'))
summary = summary[!duplicated(summary[,'projid...2']),]
rownames(summary) = summary[['projid...2']]

df = t(as.data.frame(av_expression$Oli[c('SOAT1','SOAT2','APOE'),]))
df = cbind(df, summary[as.character(rownames(df)),c('apoe_genotype','niareagansc','cogdx')])
df$AD = ifelse(df$niareagansc%in%c(1,2) & df$cogdx==4, 'AD','nonAD')
df$APOE_geno = ifelse(df$apoe_genotype==33,'APOE3/3', 'APOE3/4 & APOE4/4')

pdf('../plots/APOE_oligo_sc.pdf', width = 4, height = 4)
ggplot(df, aes(x=as.factor(APOE_geno), y=APOE)) +
 geom_boxplot() + facet_wrap(~AD, ncol = 2) + geom_jitter()+ stat_compare_means(method = "wilcox.test", comparisons = list(c('APOE3/3', 'APOE3/4 & APOE4/4')))
dev.off()
# don't plot SOAT2, as barely detected in oligosSpeicherXI282

 #ggplot(df, aes(x=as.factor(APOE), y=SOAT2)) +
  #geom_boxplot() + facet_wrap(~AD, ncol = 2) + geom_jitter()+ stat_compare_means(method = "wilcox.test", comparisons = list(c('APOE3/3', 'APOE3/4 & APOE4/4')))

usage = assays(e)$usage
# get usage
usage['SOAT1','Oli']

usage['SOAT2','Oli']

# also show SOAT1 and SOAT2 distributions for the bulk APOE3 and 4 RNA sequencing data

norm_counts_no_drug = read.csv('../data_outputs/ipsc_bulk_no_drug.csv')
x = as.data.frame(read.table('../raw_data/191210Tsa/RNAseqQC-4367S/D19-13224-4367S_geneexp.txt', header = FALSE))
gene.df = as.data.frame(cbind(as.character(x$V1), as.character(x$V3)))
rownames(gene.df) = gene.df$V1

meta1 = read.table('../raw_data/APOE_ipsc_opc_sequencing_meta.csv', sep = ',', header = F)
var = as.character(meta1$V2)
var = ifelse(var=='APOE4', 1, 0)

df = t(as.data.frame(norm_counts_no_drug[gene.df[norm_counts_no_drug$X, 'V2']=='SOAT1',!colnames(norm_counts_no_drug)%in%c('X')]))
colnames(df) = c('SOAT1')
df = as.data.frame(df)
df$APOE = var

ggplot(df, aes(x=as.factor(APOE), y=SOAT1)) +
 geom_boxplot() + geom_jitter()+ stat_compare_means(method = "wilcox.test", comparisons = list(c('0', '1')))


# SOAT2 values are all zero
# unname(norm_counts_no_drug[gene.df[norm_counts_no_drug$X, 'V2']=='SOAT2',]))))

# also look at APOE expression
df = t(as.data.frame(norm_counts_no_drug[gene.df[norm_counts_no_drug$X, 'V2']=='APOE',!colnames(norm_counts_no_drug)%in%c('X')]))
colnames(df) = c('APOE')
df = as.data.frame(df)
df$APOE_geno = var

pdf('../plots/APOE_ipsc_oligo.pdf', width = 2, height = 2)
ggplot(df, aes(x=as.factor(APOE_geno), y=APOE)) +
 geom_boxplot() + geom_jitter()+ stat_compare_means(method = "wilcox.test", comparisons = list(c('0', '1')))
dev.off()
