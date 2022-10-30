source('../functions/qc_and_annotation_aux_functions.r')
library(SingleCellExperiment)
library(tidyr)

print('reading sce')
sce = readRDS('../data/single_cell_data/single_cell_experiment_object.rds')
meta = colData(sce)
cell_labels = rownames(meta)

# in how many cells per celltype is each gene detected?
print('getting nonzero counts')
counts_nonzero = counts(sce)>0
detected_genes_cell = sum_counts(counts_nonzero, label = as.data.frame(meta$cell.type), cell_labels)
fraction_detected_genes_cell = t(apply(detected_genes_cell$summed_counts, 1, function(x){x/detected_genes_cell$ncells}))

# get expression list
fraction_detected_genes_cell_binary = fraction_detected_genes_cell>.10
genes = rownames(fraction_detected_genes_cell)
expressed = lapply(colnames(fraction_detected_genes_cell_binary), function(x) genes[unname(fraction_detected_genes_cell_binary[,x])])
names(expressed) = colnames(fraction_detected_genes_cell_binary)

print('saving expressed genes')
saveRDS(expressed, '../data/single_cell_data/expressed_genes_per_celltype.rds')
