########## effects stratification for extended data figure 3 #############
##########################################################################

print('|| stratification of E4 effect... ||')

# load functions
source('../functions/pathway_fits.r')

# required packages
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
library('ggpubr')

# load the pathways
data = readRDS('../data/other_analyses_outputs/pathways.rds')
low_removed_bp = data$pathways$low_removed_bp
apoe_gsets_low_removed = data$pathways$apoe_gsets_low_removed
apoe_gsets = data[['pathways']][['apoe_gsets_all']]
all_paths = data$pathways$all

# output data
all_data = list()

# load the data
expressed = readRDS('../data/single_cell_data/expressed_genes_per_celltype.rds')
av_expression = readRDS('../data/single_cell_data/individual_level_averages_per_celltype.rds')
summary = read.csv('../data/single_cell_data/metadata_by_individual.csv')
summary$APOE4 = ifelse(summary$apoe_genotype == '33','e3','e4')
summary$AD = ifelse(summary$niareagansc%in%c(1,2),1,0)
rownames(summary) = summary[,'projid...2']
low_removed_lipid_paths = data$pathways$low_removed_lipid_associated

# run GSVA on lipid pathways
order = c('Oli')
out_lipid_terms = lapply(order, function(x) t(gsva(as.matrix(av_expression[[x]]), low_removed_lipid_paths[[x]], mx.diff=TRUE, verbose=F, kcdf=c("Gaussian"), min.sz=5, max.sz = 150, parallel.sz=50)))
names(out_lipid_terms) = order
all_data[['res']][['lipid_associated']][['gsva_out']] = out_lipid_terms

# E4 effect, stratified by AD status
meta = summary[summary$niareagansc%in%c(3,4)& summary$apoe_genotype!=44,]
mod = model.matrix(~APOE4 + amyloid + nft + age_death + msex + pmi, data=meta)
all_data[['res']][['lipid_associated']][['APOE34_effect_nia34']] = get_stratified_fits(meta, out_lipid_terms, mod, 'APOE4e4')

meta = summary[summary$niareagansc%in%c(1,2)& summary$apoe_genotype!=44,]
mod = model.matrix(~APOE4 + amyloid + nft + age_death + msex + pmi, data=meta)
all_data[['res']][['lipid_associated']][['APOE34_effect_nia12']] = get_stratified_fits(meta, out_lipid_terms, mod, 'APOE4e4')

# AD effect, stratified by APOE4
mod = model.matrix(~APOE4 + AD + age_death + msex + pmi, data=summary)
all_data[['res']][['lipid_associated']][['AD_effect_unstratified']] = get_stratified_fits(summary, out_lipid_terms, mod, 'AD')

meta = summary[summary$apoe_genotype==33,]
mod = model.matrix(~AD + age_death + msex + pmi, data=meta)
all_data[['res']][['lipid_associated']][['AD_effect_apoe33']] = get_stratified_fits(meta, out_lipid_terms, mod, 'AD')

meta = summary[summary$apoe_genotype==34,]
mod = model.matrix(~AD + age_death + msex + pmi, data=meta)
all_data[['res']][['lipid_associated']][['AD_effect_apoe34']] = get_stratified_fits(meta, out_lipid_terms, mod, 'AD')

# get all the effect sizes and p-values as a table for pathway of interest
cholest_path = 'cholesterol biosynthesis III (via desmosterol) Homo sapiens PWY66-4'
out = list()
names = names(all_data[['res']][['lipid_associated']])
names = names[!names%in%c('gsva_out')]

for(i in names){
    f = all_data[['res']][['lipid_associated']][[i]][cholest_path,c('logFC', 'P.Value')]
    f$name = i
    out[[i]] = f
}

# save the table with logFC values and p-value
write.csv(do.call('rbind', out),'../data/supplementary_tables/e4_stratified_stats.csv')
saveRDS(all_data, '../data/other_analyses_outputs/stratified_anaylsis.rds')
print('done')
