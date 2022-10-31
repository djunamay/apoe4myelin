##############################################################################################################
library(nebula)
library(SingleCellExperiment)
print('loading sce')
sce = readRDS('../data/single_cell_data/single_cell_experiment_object_qced.rds')
print('subsetting oli')
sce.oli <- sce[,sce$cell.type%in%"Oli"]
gc()
##############################################################################################################
print('setting up model')
C <- counts(sce.oli)
Predictors <- as.data.frame(sce.oli@colData[,c("amyloid","nft", "APOE4", "age_death", "batch")])
Predictors$amyloid <- scale(Predictors$amyloid)
Predictors$nft <- scale(Predictors$nft)
Predictors$age_death <- scale(Predictors$age_death)
Predictors$APOE4 <- as.character(Predictors$APOE4)
Predictors$batch <- as.character(Predictors$batch)
df = model.matrix(~APOE4+amyloid+nft+age_death+batch, data=Predictors)
print('computing degs')
sparsematrix <- as(C, "CsparseMatrix")

re = nebula(sparsematrix, sce.oli$projid, pred=df, offset=Matrix::colSums(sparsematrix))
##############################################################################################################
out <- re$summary
rownames(out) <- out$gene
out$score <- -log10(out[[i]]$p_APOE4yes) * ifelse(out[[i]]$logFC_APOE4yes>0, 1, -1)
print('saving degs')
saveRDS(out, file="../data/differentially_expressed_genes/E4_nebula_associations_Oli.rds")
##############################################################################################################
