# load functions
source('../functions/pathway_analyses.r')
library(ComplexHeatmap)
library(GSVA)

# load the data
data = readRDS('../data/pathways.rds')
low_removed_bp = data$pathways$low_removed_bp
apoe_gsets_low_removed = data$pathways$apoe_gsets_low_removed
apoe_gsets = data[['pathways']][['apoe_gsets_all']]
all_paths = data$pathways$all

# output data
all_data = list()

set.seed(5)

# load libraries
print('loading libraries..')
library("readxl")
library('SingleCellExperiment')
library('GSVA')
library('limma')
library('parallel')
library('tidyr')
library(stringr)
library(GSA)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)

# load the data
print('loading the data..')
e = readRDS('../../submission_code_09012021/data/Summary.DE.celltype.rds')
expressed = metadata(e)$expressed.genes
av_expression = readRDS('../../submission_code_09012021/data/Averages.by.celltype.by.individual.rds')

summary = as.data.frame(read_excel('../../submission_code_09012021/data/Overview_PFC_projids_data_384.xlsx'))
summary = summary[!duplicated(summary[,'projid...2']),]
rownames(summary) = summary[['projid...2']]

# load the data from figure 1
f1_data = readRDS('../data/data_figure_1.rds')
f2_data = readRDS('../data/figure_2_data.rds')

geno = summary[colnames(df), 'apoe_genotype']

oli.lipid.fits = f1_data$res$lipid_associated$fits$Oli
lipid.paths.oli = f1_data$pathways$lipid_associated$Oli
paths = rownames(oli.lipid.fits[oli.lipid.fits$P.Value<0.05,])
biosynth_genes = unique(unname(unlist(lipid.paths.oli[paths])))

###### save the heatmap

# save plot for each of the signature
genes = c('DHCR24', 'LPIN2', 'IRS2', 'NR1H2', 'LBR', 'LPIN1')
genes1 = c('PCYT1B', 'SEC23A', 'PRKN', 'SCP2', 'LPCAT3')
genes2 = c('HSP90B1', 'CALR', 'HSPA5', 'MBTPS1', 'ATF6')
o = list()
o[['sig1']] = genes
o[['sig2']] = genes1
o[['sig3']] = genes2
o[['sig4']] = biosynth_genes

for(i in names(o)){
  df = av_expression$Oli[c(o[[i]]),]
  df = t(scale(t(df)))

  anns = ifelse(summary[colnames(df), 'apoe_genotype']==33, 'E3', 'E4')
  anns2 = ifelse(summary[colnames(df), 'niareagansc']%in%c(1,2), 'AD', 'nonAD')

  column_ha = HeatmapAnnotation(geno = anns)
  column_ha2 = HeatmapAnnotation(geno = anns2)

  ha = rowAnnotation(e3 = anno_density(df[,anns=='E3'], xlim = c(-2,4)))
  ha2 = rowAnnotation(e4 = anno_density(df[,anns=='E4'], xlim = c(-2,4)),
  gp = gpar(fill = "yellow"))

  #pdf(paste0('../plots/pseudobulk_signature_',i,'.pdf'), width=7, height = 3)
  print(Heatmap(df, bottom_annotation = c(column_ha,column_ha2), left_annotation = c(ha, ha2)))
  #dev.off()
}
