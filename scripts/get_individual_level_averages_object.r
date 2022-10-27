library(tidyr)
library(ACTIONet)
source('../functions/qc_and_annotation_aux_functions.r')
sce = readRDS('../data/single_cell_data/single_cell_experiment_object_qced.rds')
sce =  scran.normalize(sce)
meta = colData(sce)
cell_labels = rownames(meta)

# use matrix multiplication to summarize (sum) across counts per cell (including all individuals)
labels = as.data.frame(as.character(interaction(meta$cell.type, meta$projid)))
summed_logcounts_cellxind = sum_counts(logcounts(sce), labels, cell_labels)

# get averages corresponding to both count matrices
avs_logcounts_cellxind = t(apply(summed_logcounts_cellxind$summed_counts, 1, function(x){x/summed_logcounts_cellxind$ncells}))

# split the averages by celltype
x = (strsplit(colnames(avs_logcounts_cellxind), '[.]'))
celltype = unlist(lapply(1:length(x), function(i) x[[i]][[1]]))
individual = unlist(lapply(1:length(x), function(i) x[[i]][[2]]))
celltype_unique = unique(celltype)
avs_by_ind_out = list()
for(i in celltype_unique){
    index = celltype==i
    df = avs_logcounts_cellxind[, index]
    colnames(df) = individual[index]
    avs_by_ind_out[[i]] = df
}
saveRDS(avs_by_ind_out, '../data/single_cell_data/individual_level_averages_per_celltype2.rds')
