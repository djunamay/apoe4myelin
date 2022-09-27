########## pseudo-bulk visualization for extended data figure 2 ##########
##########################################################################

print('|| pseudobulk visualizations for pathways of interest... ||')

#required packages
print('loading libraries..')
library("readxl")
library('SingleCellExperiment')
library('GSVA')
library('limma')
library('parallel')
library('tidyr')
library('stringr')
library('GSA')
library('ComplexHeatmap')
library('circlize')
library('ggplot2')

# load the data
expressed = readRDS('../data/single_cell_data/expressed_genes_per_celltype.rds')
av_expression = readRDS('../data/single_cell_data/individual_level_averages_per_celltype.rds')
summary = read.csv('../data/single_cell_data/metadata_by_individual.csv')
summary$APOE4 = ifelse(summary$apoe_genotype == '33','e3','e4')
rownames(summary) = summary[,'projid']

# load the data from figure 1
cholest = readRDS('../data/other_analyses_outputs/cholesterol_analysis.rds')
geno = summary[colnames(df), 'apoe_genotype']

# save plot for each of the signature
o = list()
o[['sig1']] = cholest$union_cholest_biosynth$genes

for(i in names(o)){
  df = av_expression$Oli[c(o[[i]]),]
  df = t(scale(t(df)))

  anns = ifelse(summary[colnames(df), 'apoe_genotype']==33, 'E3', 'E4')
  column_ha = HeatmapAnnotation(APOE_genotype = anns)

  ha = rowAnnotation(e3 = anno_density(df[,anns=='E3'], xlim = c(-2,4)))
  ha2 = rowAnnotation(e4 = anno_density(df[,anns=='E4'], xlim = c(-2,4)))

  pdf(paste0('../plots/Extended_2/pseudobulk_signature_',i,'.pdf'), width=15, height = 5)
  print(Heatmap(df, bottom_annotation = c(column_ha), left_annotation = c(ha, ha2)))
  dev.off()
}

print('done.')
