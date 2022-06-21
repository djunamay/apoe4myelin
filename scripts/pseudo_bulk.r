# load functions
source('../functions/pathway_analyses.r')
library(ComplexHeatmap)
library(GSVA)

# load the data
data = readRDS('../data/other_analyses_outputs/pathways.rds')
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
expressed = readRDS('../data/single_cell_data/expressed_genes_per_celltype.rds')
av_expression = readRDS('../data/single_cell_data/individual_level_averages_per_celltype.rds')
summary = read.csv('../data/single_cell_data/metadata_by_individual.csv')
summary$APOE4 = ifelse(summary$apoe_genotype == '33','e3','e4')
rownames(summary) = summary[,'projid...2']

print('running')
# load the data from figure 1
cholest = readRDS('../data/other_analyses_outputs/cholesterol_analysis.rds')
er_stress = readRDS('../data/other_analyses_outputs/er_stress_results.rds')

geno = summary[colnames(df), 'apoe_genotype']

# save plot for each of the signature
o = list()
o[['sig1']] = cholest$union_cholest_biosynth$genes
o[['sig2']] = names(er_stress$ER$ATF6_pathway_degs)

for(i in names(o)){
  df = av_expression$Oli[c(o[[i]]),]
  df = t(scale(t(df)))

  anns = ifelse(summary[colnames(df), 'apoe_genotype']==33, 'E3', 'E4')
  anns2 = ifelse(summary[colnames(df), 'niareagansc']%in%c(1,2), 'AD', 'nonAD')

  column_ha = HeatmapAnnotation(APOE_genotype = anns)
  column_ha2 = HeatmapAnnotation(NIAreagansc = anns2)

  ha = rowAnnotation(e3 = anno_density(df[,anns=='E3'], xlim = c(-2,4)))
  ha2 = rowAnnotation(e4 = anno_density(df[,anns=='E4'], xlim = c(-2,4)))

  #pdf(paste0('../plots/pseudobulk_signature_',i,'.pdf'), width=7, height = 3)
  print(Heatmap(df, bottom_annotation = c(column_ha,column_ha2), left_annotation = c(ha, ha2)))
  #dev.off()
}

print('done.')
