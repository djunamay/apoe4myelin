########## compute wilcoxon degs ##########
###########################################
library(SingleCellExperiment)
print('loading the sce object')
source('../functions/differential_expression.r')
source('../functions/qc_and_annotation_aux_functions.r')

sce = readRDS('../data/single_cell_data/single_cell_experiment_object_qced.rds')

# for all
print('subsetting sce object')
oli = sce[,sce$cell.type=='Oli']
print('computing degs')
out_all = Wilcox.differential(logcounts(oli), ifelse(oli$apoe_genotype==33, 0, 1))

# stratified by pathology
oli_ad = oli[,oli$AD=='yes' & oli$apoe_genotype!=44]
out_ad = Wilcox.differential(logcounts(oli_ad), ifelse(oli_ad$apoe_genotype==33, 0, 1))

oli_nonad = oli[,oli$AD=='no' & oli$apoe_genotype!=44]
out_nonad = Wilcox.differential(logcounts(oli_nonad), ifelse(oli_nonad$apoe_genotype==33, 0, 1))

saveRDS(out_all, '../data/differentially_expressed_genes/oli_wilcox_results.rds')
saveRDS(out_ad, '../data/differentially_expressed_genes/oli_wilcox_results_AD.rds')
saveRDS(out_nonad, '../data/differentially_expressed_genes/oli_wilcox_results_noAD.rds')

print('done.')
