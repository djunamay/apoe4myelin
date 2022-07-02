########## compute wilcoxon degs ##########
###########################################

source('../functions/differential_expression.r')

sce = readRDS('../data/single_cell_data/single_cell_experiment_object.rds')

# for all
oli = sce[,sce$cell.type=='Oli']
out_all = Wilcox.differential(logcounts(oli), ifelse(oli$apoe_genotype==33, 0, 1))

# stratified by pathology
oli_ad = oli[,oli$niareagansc%in%c(1,2) & oli$apoe_genotype!=44]
out_ad = Wilcox.differential(logcounts(oli_ad), ifelse(oli_ad$apoe_genotype==33, 0, 1))

oli_nonad = oli[,oli$niareagansc%in%c(3,4) & oli$apoe_genotype!=44]
out_nonad = Wilcox.differential(logcounts(oli_nonad), ifelse(oli_nonad$apoe_genotype==33, 0, 1))

saveRDS(out_all, '../data/differentially_expressed_genes_data/oli_wilcox_results.rds')
saveRDS(out_ad, '../data/differentially_expressed_genes_data/oli_wilcox_results_AD.rds')
saveRDS(out_nonad, '../data/differentially_expressed_genes_data/oli_wilcox_results_noAD.rds')

print('done.')
