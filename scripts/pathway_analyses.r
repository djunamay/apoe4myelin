########## run pathway analyses (GSVA) ##########
#################################################
print('|| running GSVA on pathway databases... ||')

##### required libraries ####
source('../functions/pathway_fits.r')
library('GSVA')
library("readxl")
library('SingleCellExperiment')
library('limma')
library('parallel')
library('tidyr')
library('stringr')

# load the pathways
pathways = readRDS('../data/other_analyses_outputs/pathways.rds')
low_removed_bp = pathways$pathways$low_removed_bp
apoe_gsets_low_removed = pathways$pathways$apoe_gsets_low_removed
apoe_gsets = pathways$pathways$apoe_gsets_all
all_paths = pathways$pathways$all
low_removed_lipid_paths = pathways$pathways$low_removed_lipid_associated

# load the data
expressed = readRDS('../data/single_cell_data/expressed_genes_per_celltype.rds')
av_expression = readRDS('../data/single_cell_data/individual_level_averages_per_celltype.rds')
summary = read.csv('../data/single_cell_data/metadata_by_individual.csv')
summary$APOE4 = ifelse(summary$apoe_genotype == '33','e3','e4')
rownames(summary) = summary[,'projid...2']

all_data = list()

# perform enrichment analyses for different pathway databases
print('for GO BP pathways...')
order = c('In', 'Ex', 'Oli', 'Opc', 'Mic','Ast')
print('GO BP...')
out_bp = get_pathway_fits(order, av_expression, low_removed_bp, top_20=TRUE, summary)
all_data[['GO_BP']] = out_bp

print('for APOE-associated pathways...')
out_apoe = get_pathway_fits(order, av_expression, apoe_gsets_low_removed, top_20=FALSE, summary)
all_data[['APOE']] = out_apoe

print('for lipid-associated pathways...')
out_lipid = get_pathway_fits(order, av_expression, low_removed_lipid_paths, top_20=FALSE, summary)
all_data[['lipid']] = out_lipid

# save all the data
print('saving the data...')
saveRDS(all_data, '../data/other_analyses_outputs/pathway_scores.rds')
print('done.')
