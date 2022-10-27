##############################################################################################################
library(nebula)
sce = readRDS('../data/single_cell_data/single_cell_experiment_object.rds')
sce.oli <- sce[,sce$cell.type%in%"Oli"]
gc()
##############################################################################################################
C <- counts(sce.oli)
Predictors <- as.data.frame(sce.oli@colData[,c("amyloid","nft", "APOE4", "age_death", "batch")])
Predictors$amyloid <- scale(Predictors$amyloid)
Predictors$nft <- scale(Predictors$nft)
Predictors$age_death <- scale(Predictors$age_death)
Predictors$APOE4 <- as.character(Predictors$APOE4)
Predictors$batch <- as.character(Predictors$batch)
df = model.matrix(~APOE4+amyloid+nft+age_death+batch, data=Predictors)
re = nebula(C, sce.oli$projid, pred=df, offset=Matrix::colSums(C))
##############################################################################################################
out <- re$summary
rownames(out) <- out$gene
out$score <- -log10(DEoutList[[i]]$p_APOE4yes) * ifelse(DEoutList[[i]]$logFC_APOE4yes>0, 1, -1)
saveRDS(out, file="../APOE4_Oligo_project/data/E4_nebula_associations_Oli.rds")
##############################################################################################################
