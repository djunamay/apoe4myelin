#################################################################################
# run analysis
library(ACTIONet)
library(BiocParallel)
library(NetLibR)
library(ComplexHeatmap)
#################################################################################
source("../functions/qc_and_annotation_aux_functions.r")
#################################################################################
# load raw snRNAseq counts in single cell experiment format
sce <- readRDS("../data/single_cell_data/single_cell_experiment_object_raw.rds")
#################################################################################
# basic QC filtering 
sce <- QC.filter.sce(sce)
#[1] "initial dim"
#[1]  33538 178789
#[1] "cells depth"
#[1] 178789 175137
#[1] "genes non detected"
#[1] 1740
#[1] "cells mito"
#[1] 1151
#[1] "final dim"
#[1]  17915 173986
gc()
##############################################################################################################
#run clustering
sce <- reduce.and.batch.correct.sce.Harmony(sce, batch.vec = as.character(sce$batch))
sce <- ACTIONet.clustering.analysis(sce) # had to comment out Extract.tops.by.column in function file
##############################################################################################################
# cluster marker-based annotation
RefCellTypeMarkers <- readRDS("../data/single_cell_data/RefCellTypeMarkers.adultBrain.rds")
PanglaoDB <- readRDS("../data/single_cell_data/PanglaoDB.by.organ.by.celltype.rds")

sce <- Add.best.annotation(sce, geneSet = RefCellTypeMarkers, Annotlabel = "PsychENCODE.celltypes", ClustLabel = "cluster")
sce <- Add.best.annotation(sce, geneSet = PanglaoDB$Brain, Annotlabel = "PanglaoDB.celltypes", ClustLabel = "cluster")
##############################################################################################################
# plot diagnostic plots for manual QC rounds
Cell.colors <- c(Ast="#E64B35FF", Ex="#4DBBD5FF", In="#00A087FF", Mic="#3C5488FF", Oli="#F39B7FFF", Opc="#8491B4FF", Endo="#FF7F00", Per="#BC80BD", na="lightgrey", macroglia="#C4451C", immune="#683B79", neuronal="#66B0FF", vascular="#7E6148FF", other="#E4E1E3", Tcell="#DEA0FD",Macrophage="#683B79",Fib="#FEAF16")

pdf("../plots/qc_annotation/cluster.stats.pdf", width = 6, height = 10)
Plot.cluster.stats.annot(sce, cell.type.colors = Cell.colors, varNames = "cluster") # only works with  some of the elements commented out in the function file
dev.off()
Openfile("cluster.stats.pdf")

pdf("../plots/qc_annotation/cluster.top.genes.pdf", height = 10)
Plot.top.genes(sce)
dev.off()
Openfile("../plots/qc_annotation/cluster.top.genes.pdf")
##############################################################################################################
sce.out <- Extract.summarized.pattern.output.from.sce(sce, getEnrichments = T)

Map <- sce.out$enrichments$PsychENCODE.celltypes.enrichment$annotation[order(sce.out$enrichments$PsychENCODE.celltypes.enrichment$annotation$labels),]
Corder <- rownames(Map)
Annot <- sce.out$enrichments$PsychENCODE.celltypes.enrichment$annotation[Corder, "labels"]
Annot <- rowAnnotation(cell.group=Annot, col=list(cell.group=Cell.colors))
Annot2 <- rowAnnotation(count=anno_barplot(as.numeric(table(sce.out$colMetadata$cluster)[Corder])))

pdf("../plots/qc_annotation/pre_celltype_enrichments.pdf", width = 10, height = 6)
Heatmap(Filter.significant.enrichment.map.values(sce.out$enrichments$PsychENCODE.celltypes.enrichment$enrichment)[rownames(Map),unique(Map$labels)], rect_gp = gpar(col = "black"), row_names_side = "left", name="z-score\n(PsychEncode)", cluster_rows = F, cluster_columns = F) +
  Heatmap(Filter.significant.enrichment.map.values(sce.out$enrichments$PanglaoDB.celltypes.enrichment$enrichment[rownames(Map),]), rect_gp = gpar(col = "black"), name="z-score\n(other)") + 
  Annot +
  Annot2
dev.off()
Openfile("../plots/qc_annotation/pre_celltype_enrichments.pdf")
##############################################################################################################
pdf("../plots/qc_annotation/pre_2D_cells_clusters.pdf")
Plot.coords(as.data.frame(sce@colData), Col = Vector.to.colors(sce@colData$cluster), add.text = T, add.legend = F)
dev.off()
Openfile("../plots/qc_annotation/pre_2D_cells_clusters.pdf")

pdf("../plots/qc_annotation/pre_2D_cells_groups.pdf")
Plot.coords(as.data.frame(sce@colData), Col = Cell.colors[sce@colData$PsychENCODE.celltypes], add.text = T, add.legend = F, Trans = sce@colData$connectivity)
dev.off()
Openfile("../plots/qc_annotation/pre_2D_cells_groups.pdf")

pdf("../plots/qc_annotation/pre_2D_cells_celltypes.pdf") # not working
Plot.coords(as.data.frame(sce@colData), Col = Cell.colors[sce@colData$PsychENCODE.celltypes], add.text = T, add.legend = F, Trans = sce@colData$connectivity)
dev.off()
Openfile("../plots/qc_annotation/pre_2D_cells_celltypes.pdf")

pdf("../plots/qc_annotation/pre_2D_cells_celltypes.PanglaoDB.pdf") # not working
Plot.coords(as.data.frame(sce@colData), Col = Vector.to.colors(sce@colData$PanglaoDB.celltypes), add.text = T, add.legend = F, Trans = sce@colData$connectivity)
dev.off()
Openfile("../plots/qc_annotation/pre_2D_cells_celltypes.PanglaoDB.pdf")
##############################################################################################################
toRemove <- c("C14")
sum(sce$cluster%in%toRemove) # 3271 cells 
which(!sce$cluster%in%toRemove) -> Passed
#length(Passed) 170715
#dim(sce)
#[1]  17915 173986
sce <- Sub.sce.cols(InSCE = sce, xcIndx = Passed)
gc()
##############################################################################################################
sce <- Add.subclustering.analysis(xsce = sce, groupVar = "PsychENCODE.celltypes", BatchCorrect = TRUE) # not working
##############################################################################################################
pdf("../plots/qc_annotation/pre_subcluster.correlation.pdf", width = 13, height = 11)
Plot.subcluster.correlation(xsce = sce, xColors = Cell.colors) # could not find function "Extract.subcluster.global.signatures"
dev.off()
Openfile("../plots/qc_annotation/pre_subcluster.correlation.pdf")

pdf("../plots/qc_annotation/pre_subcluster.stats.pdf", width = 12, height = 10)
Plot.cluster.stats.annot(sce, cell.type.colors = Cell.colors, clusterLab = "sub.cluster", varNames = c("sub.cluster"))
dev.off()
Openfile("../plots/qc_annotation/pre_subcluster.stats.pdf")

pdf("../plots/qc_annotation/pre_subcluster.networks.pdf") # not working
Plot.subcluster.networks(xsce = sce, cex=0.5)
dev.off()
Openfile("../plots/qc_annotation/pre_subcluster.networks.pdf")

pdf("../plots/qc_annotation/pre_subcluster.top.genes.pdf", height = 8)
Plot.top.genes.all.subclusters(xsce = sce)
dev.off()
Openfile("../plots/qc_annotation/pre_subcluster.top.genes.pdf")

pdf("../plots/qc_annotation/pre_subcluster.phylo.pdf", height = 12)
Plot.subcluster.phylo(xsce = sce, addColor = TRUE) # could not find function "Extract.subcluster.global.signatures"
dev.off()
Openfile("../plots/qc_annotation/pre_subcluster.phylo.pdf")
##############################################################################################################
saveRDS(sce, file = "../data/single_cell_data/sce.rds")
##############################################################################################################
pdf("../plots/qc_annotation/pre_2D_cells_clusters.pdf")
Plot.coords(sce@colData, Col = Vector.to.colors(sce@colData$cluster), add.text = T, add.legend = F)
dev.off()
Openfile("../plots/qc_annotation/pre_2D_cells_clusters.pdf")

pdf("../plots/qc_annotation/pre_2D_cells_groups.pdf") # not working
Plot.coords(sce@colData, Col = Cell.colors[sce@colData$PsychENCODE.celltypes], add.text = T, add.legend = F, Trans = sce@colData$connectivity)
dev.off()
Openfile("../plots/qc_annotation/pre_2D_cells_groups.pdf")

pdf("../plots/qc_annotation/pre_2D_cells_celltypes.pdf") # not working
Plot.coords(sce@colData, Col = Cell.colors[sce@colData$PsychENCODE.celltypes], add.text = T, add.legend = F, Trans = sce@colData$connectivity)
dev.off()
Openfile("../plots/qc_annotation/pre_2D_cells_celltypes.pdf")

pdf("../plots/qc_annotation/pre_2D_cells_celltypes.PanglaoDB.pdf") # not working
Plot.coords(sce@colData, Col = Vector.to.colors(sce@colData$PanglaoDB.celltypes), add.text = T, add.legend = F, Trans = sce@colData$connectivity)
dev.off()
Openfile("../plots/qc_annotation/pre_2D_cells_celltypes.PanglaoDB.pdf")
##############################################################################################################
# visualize putative doublet and low quality cells in 2D plot
toRemove <- c("Ast.C6", "Ast.C7", "Ast.C5", "Ex.C14", "Ex.C9", "In.C12", "Oli.C6", "Opc.C5", "Opc.C6", "Opc.C7", "Opc.C8", "Opc.C9", "Mic.C7", "Mic.C8", "Mic.C9", "Endo.C5", "Endo.C6")
toKeep <- which(!sce@colData$sub.cluster%in%toRemove)  # not working
Plot.coords(sce@colData) # not working
Plot.coords(sce@colData, Col = ifelse(sce@colData$sub.cluster%in%toRemove, "red", "grey")) # not working
Plot.coords(sce@colData[which(sce@colData$sub.cluster%in%toRemove),]) # not working
Plot.coords(sce@colData[toKeep,]) # not working

# remove subclustering + outlier cells
Oulier.cells <- Get.outlier.cells(sce@colData[toKeep,])
toKeep <- which(!sce@colData$sub.cluster%in%toRemove)
toKeep <- which( (!sce@colData$sub.cluster%in%toRemove) & (!colnames(sce)%in%Oulier.cells) )
# visualize filtered cells
Plot.coords(sce@colData)
Plot.coords(sce@colData[toKeep,])
##############################################################################################################
sce <- Sub.sce.cols(sce, xcIndx = toKeep)
dim(sce)
#dim: [1]  17915 164741
gc()
##############################################################################################################
xx <- Get.reLayout.corrdinates(sce) # not working
sce@colData <- DataFrame(sce@colData, xx)
Plot.coords(xx, x = "re.x", y = "re.y", Col = Cell.colors[sce$PsychENCODE.celltypes], Trans = sce$connectivity)
sce@metadata$subclustering.out <- NULL
sce@metadata <- list()
gc()

saveRDS(sce, file="../data/single_cell_data/sce2.rds")
##############################################################################################################
Get.top.matching.PanglaoDB(sce@metadata$subclustering.out$Endo@metadata$cluster.signatures) # cant find this function
#Mic.C6 - T cells 
#Endo.C1 - Endo
#Endo.C2 - Per
#Endo.C3 - Fibro
#Endo.C4 - SMC
##############################################################################################################
Mapping <- as.data.frame(unique(sce@colData[,c("PsychENCODE.celltypes", "sub.cluster")]))
rownames(Mapping) <- Mapping$sub.cluster # not working (non-unique values)
Mapping$celltype.relable <- Mapping$PsychENCODE.celltypes
Mapping$celltype.relable[Mapping$sub.cluster=="Mic.C6"] <- "Tcell"
Mapping$celltype.relable[Mapping$sub.cluster=="Endo.C1"] <- "Endo"
Mapping$celltype.relable[Mapping$sub.cluster=="Endo.C2"] <- "Per"
Mapping$celltype.relable[Mapping$sub.cluster=="Endo.C3"] <- "Fib"
Mapping$celltype.relable[Mapping$sub.cluster=="Endo.C4"] <- "SMC"
View(Mapping)
Mapping$celltype.colors <- All.colors[Mapping$celltype.relable]
saveRDS(Mapping, file="../data/single_cell_data/Mapping.rds")
##############################################################################################################
# update cell labels
Mapping <- readRDS("../data/single_cell_data/Mapping.rds")
mean(unique(sce$sub.cluster)%in%rownames(Mapping))
sce$cell.type <- Mapping[sce$sub.cluster,"celltype.relable"]
##############################################################################################################
sce@colData$celltype.relable <- Mapping[sce$sub.cluster,"celltype.relable"]
Summary.data.celltype <- Get.pseudobulk.summary(xsce = sce, groupVect = sce$celltype.relable, doParallel=TRUE) # not working
saveRDS(Summary.data.celltype, file = "../data/single_cell_data/Summary.data.celltype.rds")
saveRDS(sce, file="../data/single_cell_data/Results_Djuna/snRNA.data.sce.rds")
##############################################################################################################
