av_expression = readRDS('../data/single_cell_data/individual_level_averages_per_celltype.rds') # check are these the average logcounts? does this make sense for the
# distribution that limma expects?

source('../functions/differential_expression.r')

### required libraries
library("edgeR")
library("limma")
library('ggplot2')

summary = read.csv('../data/single_cell_data/metadata_by_individual.csv')
summary$APOE4 = ifelse(summary$apoe_genotype == '33','e3','e4')
rownames(summary) = summary[,'projid...2']
expressed = readRDS('../data/single_cell_data/expressed_genes_per_celltype.rds')


all_out = list()

for(i in names(av_expression)){
  out = RunDiffExprAnalysisLimma_pseudobulk(av_expression[[i]], 'APOE4e4', summary[colnames(av_expression[[i]]),])
  out = out$res[expressed[[i]],]
  out = out[order(out$P.Value,decreasing = F),]
  out$celltype = i
  out$gene = rownames(out)
  all_out[[i]] = out
}

pseudo_bulk = do.call('rbind', all_out)

write.csv(pseudo_bulk, '../data/other_analyses_outputs/pseudo_bulk_degs_single_cell_all_celltypes.csv')

print('done')
