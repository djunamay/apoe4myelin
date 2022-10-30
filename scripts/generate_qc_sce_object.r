##############################################################################################################
# load libraries
source('../functions/qc_and_annotation_aux_functions.r')
require(SingleCellExperiment)
require(Matrix)
##############################################################################################################
# load qc files and generate sce
print('loading matrix')
C <- readMM("../data/single_cell_data/qc_counts_data/qc_counts.mtx")
print('generating sce object')
sce_qc <- SingleCellExperiment(assays=list(counts=C))
rownames(sce_qc) <- readLines("../data/single_cell_data/qc_counts_data/qc_gene_names.txt")
colMeta <- DataFrame(read.csv("../data/single_cell_data/qc_counts_data/qc_column_metadata.csv", row.names = 1))
colnames(sce_qc) <- rownames(colMeta)
sce_qc@colData <- colMeta
print('normalizing counts')
sce_qc =  normalize.default(sce_qc)
##############################################################################################################
# plot cells 2D
plot(sce_qc@colData$x, sce_qc@colData$y, col=scales::alpha(sce_qc@colData$cell.type.color, sce_qc@colData$connectivity), cex=0.15, axes=F, xlab="", ylab="")
##############################################################################################################
print('saving')
saveRDS(sce_qc, '../data/single_cell_data/single_cell_experiment_object_qced.rds')
print('done')
