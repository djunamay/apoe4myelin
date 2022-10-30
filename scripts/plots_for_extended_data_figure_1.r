########## plots for extended data figure 1 ##########
######################################################
source('../functions/qc_and_annotation_aux_functions.r')

print('loading sce')
sce = readRDS('../data/single_cell_data/single_cell_experiment_object_qced.rds')
Mapping = readRDS('../data/single_cell_data/Mapping.rds')
Metadata = readRDS('../data/single_cell_data/Metadata.APOE.project.rds')
All.colors = readRDS('../data/single_cell_data/Cell_group_colors.rds')
Metadata$AD = as.character(Metadata$AD)

library(ComplexHeatmap)

#actionet 2D
print('actionet 2d')
pdf("../plots/Extended_1/2D_cells_celltypes_noLabs.pdf")
Cols <- Mapping[sce$sub.cluster,"celltype.colors"]
names(Cols) <- Mapping[sce$sub.cluster,"celltype.relable"]
Plot.coords(as.data.frame(sce@colData), x = "x", y = "y", Trans = sce@colData$connectivity, Col = Cols)
dev.off()

#path variables
print('path vars')
require(ggpubr)
pdf("../plots/Extended_1/Supp_PathoVariables_by_genotype_AD.pdf", height = 6)
p1 <- ggboxplot(Metadata, x = "AD", y = "amyloid") + facet_wrap(~apoe_genotype) + stat_compare_means()
p2 <- ggboxplot(Metadata, x = "AD", y = "nft") + facet_wrap(~apoe_genotype) + stat_compare_means()
ggarrange(p1, p2, nrow = 2)
dev.off()

#path variables 2
pdf("../plots/Extended_1/Supp_PathoVariables_2_by_genotype_AD.pdf", height = 6)
p1 <- ggboxplot(Metadata, x = "AD", y = "pmi") + facet_wrap(~apoe_genotype) + stat_compare_means()
p2 <- ggboxplot(Metadata, x = "AD", y = "age_death") + facet_wrap(~apoe_genotype) + stat_compare_means()
ggarrange(p1, p2, nrow = 2)
dev.off()

#marker enrichment
print('marker enrichment')
Summary.data.celltype <- readRDS("../data/single_cell_data/Summary.data.celltype.rds")
RefCellTypeMarkers <- readRDS("../data/single_cell_data/RefCellTypeMarkers.adultBrain.rds")
PanglaoDB <- readRDS("../data/single_cell_data/PanglaoDB.by.organ.by.celltype.rds")

# PsychENCODE enrichment
is.sparseMatrix <- function(x) is(x, 'sparseMatrix')
S <- Aggregate.pairwise.FC.colPairwise(assays(Summary.data.celltype)[["E"]])
x <- Permutation.enrichment.analysis(x = t(S), marker.genes = RefCellTypeMarkers)

# PanglaoDB enrichment
gSet <- Append.to.list(PanglaoDB$`Connective tissue`["Fibroblasts"], PanglaoDB$Vasculature[c(1,4)])
gSet <- Append.to.list(PanglaoDB$Brain[c("Astrocytes", "GABAergic neurons", "Glutaminergic neurons", "Microglia", "Neurons","Oligodendrocyte progenitor cells","Oligodendrocytes")], gSet)
gSet <- Append.to.list(PanglaoDB$`Immune system`["T cells"], gSet)
gSet <- Append.to.list(PanglaoDB$`Smooth muscle`["Smooth muscle cells"], gSet)
xx <- Permutation.enrichment.analysis(x = t(S), marker.genes = gSet)

temp <- x$enrichment
temp[Bonferronit(Remove.zero.pvalues(convert.z.score(temp)))>0.01] <- 0
temp <- temp[c("Ex","In","Ast","Oli","Opc","Mic","Endo","Per","Fib","SMC","Tcell"),c("Ex","In","Ast","Oli","Opc","Mic","Endo","Per")]
temp2 <- xx$enrichment
temp2[Bonferronit(Remove.zero.pvalues(convert.z.score(temp2)))>0.01] <- 0
temp2 <- temp2[rownames(temp),c("Neurons","Glutaminergic neurons","GABAergic neurons","Astrocytes","Oligodendrocytes","Oligodendrocyte progenitor cells","Microglia","Endothelial cells","Pericytes","Fibroblasts","Smooth muscle cells","T cells")]

pdf("../plots/Extended_1/Supp_Celltype_markers_enrichment.pdf", width = 5, height = 5)
Heatmap(temp, cluster_rows = F, cluster_columns = F,rect_gp = gpar(col = "black"), name="z-score a") + Heatmap(temp2, cluster_columns = F, rect_gp = gpar(col = "black"), name="z-score")
dev.off()

#cross-individual celltype correlation
print('correlations')
L = readRDS('../data/single_cell_data/individual_level_averages_per_celltype.rds')

C <- lapply(L[c("Ex", "In","Ast", "Oli", "Opc", "Mic")], function(i) cor(i)[upper.tri(cor(i))])
mean(sapply(C, mean))
#0.9435248
C <- lapply(L[c("Endo","Per","SMC","Fib","Tcell")], function(i) cor(i)[upper.tri(cor(i))])
mean(sapply(C, mean))
#0.6519454

pdf("../plots/Extended_1/Crossindividual_celltype_cor.pdf", height = 4, width = 4)
C <- lapply(L[c("Ex", "In","Ast", "Oli", "Opc", "Mic","Endo","Per","SMC","Fib","Tcell")], function(i) cor(i)[upper.tri(cor(i))])
boxplot(C, ylim=c(0,1), col=All.colors[c("Ex", "In","Ast", "Oli", "Opc", "Mic","Endo","Per","SMC","Fib","Tcell")], las=2, pch=20, cex=0.5)
abline(h=0.9435248, lty=2)
abline(h=0.6519454, lty=2)
dev.off()

#individual correlations
Lab <- rep(names(L), sapply(L, ncol))
M <- do.call(cbind, L)
C <- cor(M)
Ann <- HeatmapAnnotation(celltype=Lab, col=list(celltype=All.colors))

pdf("../plots/Extended_1/celltype_invidividual_correlation.pdf", height = 6)
Heatmap(C, bottom_annotation = Ann, show_row_names = F, show_column_names = F, name = "cor")
dev.off()

#fraction plots
print('fraction plots')
pdf("../plots/Extended_1/Supp_fraction_plots2.pdf", height = 3.5, width = 3.5)
x <- sort(apply(Table.to.matrix(table(sce$projid, sce$cell.type)), 2, median), decreasing = T)
barplot(x, col=All.colors[names(x)], las=2, ylab="median number of cells per subject")
Ns <- names(x)

#subject fraction plots
x <- sort(colMeans(Table.to.matrix(table(sce$projid, sce$cell.type))==0), decreasing = T)
barplot(x[Ns], col=All.colors[Ns], las=2, ylab="subject fraction having cells of each type", ylim=c(0,0.5))
dev.off()

#individual fractions per celltype
pdf("../plots/Extended_1/individual_fractions_perCelltype_bars.pdf", width = 5, height = 5)
barplot(Counts.to.fractions(Table.to.matrix(table(sce@colData$projid, sce@colData$cell.type))), las=2, ylab="individuals")
dev.off()

#individual fractions per celltype
dat <- Counts.to.fractions(Table.to.matrix(table(sce@colData$cell.type, sce@colData$projid)))
pdf("../plots/Extended_1/individual_fractions_perCelltype_bars_2.pdf", width = 8, height = 4)
barplot(dat, las=2, col=All.colors[rownames(dat)], ylab="individuals", xlab="")
dev.off()
