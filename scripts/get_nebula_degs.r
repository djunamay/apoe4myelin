########## plots for Figure 2 ##########
########################################
print('|| getting DEGs by NBMM... ||')

# required packages
library(nebula)

sce = readRDS('../data/single_cell_data/single_cell_experiment_object.rds')

# for all
print('computing degs...')
ace = sce[,sce$cell.type=='Oli']
ace$APOE4 = ifelse(ace$apoe_genotype%in%c('34', '44', '24'), 1, 0)

model = ~ amyloid + nft + msex + age_death + pmi + APOE4
counts = counts(ace)
meta = colData(ace)
mod = model.matrix(model, data=meta)

# run nebula
re = nebula(counts,meta$projid, pred=mod, offset=meta$NonMitoCountsg, cpc=0.1) #model = "NBGMM", kappa = 800, cutoff_cell = 20)

# remove genes that dont pass the QC

# compare nebula to the pseudo bulk results and compare nebula to the wilcox results
out = re$summary[,c('gene','logFC_APOE4', 'p_APOE4')]
out$padj = p.adjust(out$p_APOE4, method = 'fdr')
out = out[order(out$p_APOE4,decreasing = F),]

saveRDS(out, '../data/differentially_expressed_genes_data/nebula_oli_degs.rds')

print('done.')

# update with Jose's code
