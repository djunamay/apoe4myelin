########## plots for extended data figure 1 ##########
######################################################
source('../functions/qc_and_annotation_aux_functions.r')
sce = readRDS('../data/single_cell_data/single_cell_experiment_object.rds')
Mapping = readRDS('../data/single_cell_data/Mapping.rds')
Metadata = readRDS('../data/single_cell_data/Metadata.APOE.project.rds')
All.colors = readRDS('../data/single_cell_data/Cell_group_colors.rds')

#actionet 2D
pdf("../plots/Extended_1/2D_cells_celltypes_noLabs.pdf")
Cols <- Mapping[sce$sub.cluster,"celltype.colors"]
names(Cols) <- Mapping[sce$sub.cluster,"celltype.relable"]
Plot.coords(as.data.frame(sce@colData), x = "re.x", y = "re.y", Trans = sce@colData$connectivity, Col = Cols)
dev.off()

#path variables
require(ggpubr)
pdf("../plots/Extended_1/Supp_PathoVariables_by_genotype_AD.pdf", height = 6)
p1 <- ggboxplot(Metadata, x = "AD", y = "amyloid") + facet_wrap(~apoe_genotype) + stat_compare_means(comparisons = )
p2 <- ggboxplot(Metadata, x = "AD", y = "nft") + facet_wrap(~apoe_genotype) + stat_compare_means(comparisons = )
ggarrange(p1, p2, nrow = 2)
dev.off()

#path variables 2
pdf("../plots/Extended_1/Supp_PathoVariables_2_by_genotype_AD.pdf", height = 6)
p1 <- ggboxplot(Metadata, x = "AD", y = "pmi") + facet_wrap(~apoe_genotype) + stat_compare_means(comparisons = )
p2 <- ggboxplot(Metadata, x = "AD", y = "age_death") + facet_wrap(~apoe_genotype) + stat_compare_means(comparisons = )
ggarrange(p1, p2, nrow = 2)
dev.off()

#Mic, Tcell, and vascular
x <- split(sce@colData, sce@colData$PsychENCODE.celltypes)

pdf("../plots/Extended_1/Mic_Tcell_2D.pdf", width = 4, height = 4)
Plot.coords(as.data.frame(x$Mic), x = "sub.x", y = "sub.y", Col = All.colors[x$Mic$cell.type], Trans = x$Mic$sub.connectivity/max(x$Mic$sub.connectivity), cex = 0.25)
dev.off()

pdf("../plots/Extended_1/Vascular_2D.pdf", width = 4, height = 4)
Plot.coords(as.data.frame(x$Endo), x = "sub.x", y = "sub.y", Col = All.colors[x$Endo$cell.type], Trans = x$Endo$sub.connectivity/max(x$Endo$sub.connectivity), cex = 0.25)
dev.off()

#marker plots
xGenes <- c("SYT1", "NRGN", "GAD1", "VCAN", "AQP4", "MBP", "CSF1R", "CD247", "FLT1", "PDGFRB", "TAGLN", "ABCA9")
Net <- sce@metadata$network
ExpM <- assays(sce)[["logcounts"]]
Celltype.markers.imputed <- do.call(cbind, lapply(xGenes, function(i) igraph::page.rank(Net, personalized = ExpM[i,])$vector))
colnames(Celltype.markers.imputed) <- xGenes
rownames(Celltype.markers.imputed) <- colnames(ExpM)
x <- Celltype.markers.imputed
intersect(rownames(x), rownames(sce@colData)) -> keep
x <- x[keep,]
X <- sce@colData[keep,c("re.x","re.y")]

pdf("../plots/Extended_1/2D_cells_marker_plots_imputed.pdf")
for(i in colnames(x)) {
  plot(X[,"re.x"], X[,"re.y"], pch = 20, cex = 0.15, axes = F, xlab = "", ylab = "", col=scales::col_numeric(palette = blues9, domain = NULL)(x[,i]), main=i)
}
dev.off()

#marker enrichment
Summary.data.celltype <- readRDS("../data/single_cell_data/Summary.data.celltype.rds")
RefCellTypeMarkers <- readRDS("../data/single_cell_data/RefCellTypeMarkers.adultBrain.rds")
PanglaoDB <- readRDS("../data/single_cell_data/PanglaoDB.by.organ.by.celltype.rds")

# PsychENCODE enrichment
S <- Aggregate.pairwise.FC.colPairwise(Summary.data.celltype@assays[["E"]])
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

#individual correlations
L <- do.call(cbind, L)
Lab <- rep(names(L), sapply(L, ncol))
M <- do.call(cbind, L)
C <- cor(M)
Ann <- HeatmapAnnotation(celltype=Lab, col=list(celltype=All.colors))

pdf("../plots/Extended_1/celltype_invidividual_correlation.pdf", height = 6)
Heatmap(C, bottom_annotation = Ann, show_row_names = F, show_column_names = F, name = "cor")
dev.off()

#cross-individual celltype correlation
sce.by.celltype <- Split.sce.cols(sce, "cell.type")
L <- lapply(sce.by.celltype, function(i) Get.column.group.average(assays(i)[[2]], i$projid))
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

#fraction plots
All.colors <- Cell.colors
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
barplot(Counts.to.fractions(Table.to.matrix(table(sce@colData$projid, sce@colData$cell.type))), las=2, col=ACTIONet.color.bank1, ylab="individuals")
dev.off()

#individual fractions per celltype
dat <- Counts.to.fractions(Table.to.matrix(table(sce@colData$cell.type, sce@colData$projid)))
pdf("../plots/Extended_1/individual_fractions_perCelltype_bars_2.pdf", width = 8, height = 4)
barplot(dat, las=2, col=Cell.colors[rownames(dat)], ylab="individuals", xlab="")
dev.off()
