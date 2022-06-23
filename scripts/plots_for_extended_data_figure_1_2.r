##############################################################################################################
#Fig 1b
pdf("2D_cells_celltypes_noLabs.pdf")
Cols <- Mapping[sce$sub.cluster,"celltype.colors"]
names(Cols) <- Mapping[sce$sub.cluster,"celltype.relable"]
Plot.coords(as.data.frame(sce@colData), x = "re.x", y = "re.y", Trans = sce@colData$connectivity, Col = Cols)
dev.off()
Openfile("2D_cells_celltypes_noLabs.pdf")
##############################################################################################################

##############################################################################################################
#Extended data Figure 1a
require(ggpubr)
pdf("Supp_PathoVariables_by_genotype_AD.pdf", height = 6)
p1 <- ggboxplot(Metadata, x = "AD", y = "amyloid") + facet_wrap(~apoe_genotype) + stat_compare_means(comparisons = )
p2 <- ggboxplot(Metadata, x = "AD", y = "nft") + facet_wrap(~apoe_genotype) + stat_compare_means(comparisons = )
ggarrange(p1, p2, nrow = 2)
dev.off()
Openfile("Supp_PathoVariables_by_genotype_AD.pdf")
##############################################################################################################
#Extended data Figure 1b
pdf("Supp_PathoVariables_2_by_genotype_AD.pdf", height = 6)
p1 <- ggboxplot(Metadata, x = "AD", y = "pmi") + facet_wrap(~apoe_genotype) + stat_compare_means(comparisons = )
p2 <- ggboxplot(Metadata, x = "AD", y = "age_death") + facet_wrap(~apoe_genotype) + stat_compare_means(comparisons = )
ggarrange(p1, p2, nrow = 2)
dev.off()
Openfile("Supp_PathoVariables_2_by_genotype_AD.pdf")
##############################################################################################################
#Extended data Figure 1d
x <- split(sce@colData, sce@colData$PsychENCODE.celltypes)

pdf("Mic_Tcell_2D.pdf", width = 4, height = 4)
Plot.coords(as.data.frame(x$Mic), x = "sub.x", y = "sub.y", Col = All.colors[x$Mic$cell.type], Trans = x$Mic$sub.connectivity/max(x$Mic$sub.connectivity), cex = 0.25)
dev.off()

pdf("Vascular_2D.pdf", width = 4, height = 4)
Plot.coords(as.data.frame(x$Endo), x = "sub.x", y = "sub.y", Col = All.colors[x$Endo$cell.type], Trans = x$Endo$sub.connectivity/max(x$Endo$sub.connectivity), cex = 0.25)
dev.off()
##############################################################################################################
#Extended data Figure 1e
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

pdf("2D_cells_marker_plots_imputed.pdf")
for(i in colnames(x)) {
  plot(X[,"re.x"], X[,"re.y"], pch = 20, cex = 0.15, axes = F, xlab = "", ylab = "", col=scales::col_numeric(palette = blues9, domain = NULL)(x[,i]), main=i)
}
dev.off()
##############################################################################################################
#Extended data Figure 1f
Summary.data.celltype <- readRDS("Summary.data.celltype.rds")
RefCellTypeMarkers <- readRDS("RefCellTypeMarkers.adultBrain.rds")
PanglaoDB <- readRDS("PanglaoDB.by.organ.by.celltype.rds")    

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

pdf("Supp_Celltype_markers_enrichment.pdf", width = 5, height = 5)
Heatmap(temp, cluster_rows = F, cluster_columns = F,rect_gp = gpar(col = "black"), name="z-score a") + Heatmap(temp2, cluster_columns = F, rect_gp = gpar(col = "black"), name="z-score")
dev.off()
##############################################################################################################
#Extended data Figure 1g
L <- do.call(cbind, L)
Lab <- rep(names(L), sapply(L, ncol))
M <- do.call(cbind, L)
C <- cor(M)
Ann <- HeatmapAnnotation(celltype=Lab, col=list(celltype=All.colors))

pdf("celltype_invidividual_correlation.pdf", height = 6)
Heatmap(C, bottom_annotation = Ann, show_row_names = F, show_column_names = F, name = "cor")
dev.off()
Openfile("celltype_invidividual_correlation.pdf")
##############################################################################################################
#Extended data Figure 1h
sce.by.celltype <- Split.sce.cols(sce, "cell.type")
L <- lapply(sce.by.celltype, function(i) Get.column.group.average(assays(i)[[2]], i$projid))
C <- lapply(L[c("Ex", "In","Ast", "Oli", "Opc", "Mic")], function(i) cor(i)[upper.tri(cor(i))])
mean(sapply(C, mean))
#0.9435248
C <- lapply(L[c("Endo","Per","SMC","Fib","Tcell")], function(i) cor(i)[upper.tri(cor(i))])
mean(sapply(C, mean))
#0.6519454

pdf("Crossindividual_celltype_cor.pdf", height = 4, width = 4)
C <- lapply(L[c("Ex", "In","Ast", "Oli", "Opc", "Mic","Endo","Per","SMC","Fib","Tcell")], function(i) cor(i)[upper.tri(cor(i))])
boxplot(C, ylim=c(0,1), col=All.colors[c("Ex", "In","Ast", "Oli", "Opc", "Mic","Endo","Per","SMC","Fib","Tcell")], las=2, pch=20, cex=0.5)
#abline(h=mean(sapply(C, mean)), lty=2)
abline(h=0.9435248, lty=2)
abline(h=0.6519454, lty=2)
dev.off()
##############################################################################################################
#Extended data Figure 1i
All.colors <- Cell.colors
pdf("Supp_fraction_plots2.pdf", height = 3.5, width = 3.5)
x <- sort(apply(Table.to.matrix(table(sce$projid, sce$cell.type)), 2, median), decreasing = T)
barplot(x, col=All.colors[names(x)], las=2, ylab="median number of cells per subject")
Ns <- names(x)

#Extended data Figure 1j
x <- sort(colMeans(Table.to.matrix(table(sce$projid, sce$cell.type))==0), decreasing = T)
barplot(x[Ns], col=All.colors[Ns], las=2, ylab="subject fraction having cells of each type", ylim=c(0,0.5))
dev.off()

#Extended data Figure 1k
pdf("individual_fractions_perCelltype_bars.pdf", width = 5, height = 5)
barplot(Counts.to.fractions(Table.to.matrix(table(sce@colData$projid, sce@colData$cell.type))), las=2, col=ACTIONet.color.bank1, ylab="individuals")
dev.off()

#Extended data Figure 1l
dat <- Counts.to.fractions(Table.to.matrix(table(sce@colData$cell.type, sce@colData$projid)))
pdf("individual_fractions_perCelltype_bars_2.pdf", width = 8, height = 4)
barplot(dat, las=2, col=Cell.colors[rownames(dat)], ylab="individuals", xlab="")
dev.off()
##############################################################################################################


##############################################################################################################
#Extended data Figure 2c
Summary.DE.celltype <- readRDS("./share_data/Summary.DE.celltype.rds")
Celltype.averages <- assays(Summary.DE.celltype)[["E"]]

pdf("Evaluate.celltype.scores.pdf", height = 2, width = 3)
Celltype.activity <- gsva(Celltype.averages, gset.idx.list = RefCellTypeMarkers)
Celltype.activity <- Celltype.activity[colnames(Celltype.activity),]
Heatmap(Celltype.activity, cluster_rows = F, cluster_columns = F, rect_gp = gpar(col = "black"), name="geneset\nactivity score")
Map <- Table.to.matrix(table(colnames(Celltype.activity)[apply(Celltype.activity, 1, which.max)], rownames(Celltype.activity)))
Heatmap(Map, col=c("white","black"), rect_gp = gpar(col = "black"), name="assignment", cluster_rows = F, cluster_columns = F)
dev.off()
Openfile("Evaluate.celltype.scores.pdf")

pdf("Evaluate.celltype.scores.random.pdf", height = 2, width = 3)
Celltype.activity.random <- gsva(Celltype.averages, gset.idx.list = lapply(RefCellTypeMarkers, function(i) sample(rownames(sce), length(i))))
Celltype.activity.random <- Celltype.activity.random[colnames(Celltype.activity),]
Heatmap(Celltype.activity.random[colnames(Celltype.activity), rownames(Celltype.activity)], rect_gp = gpar(col = "black"), cluster_rows = F, cluster_columns = F, name="geneset\nactivity score")
Maprandom <- Table.to.matrix(table(colnames(Celltype.activity.random)[apply(Celltype.activity.random, 1, which.max)], rownames(Celltype.activity.random)))
Heatmap(Maprandom, col=c("white","black"), rect_gp = gpar(col = "black"), name="assignment", cluster_rows = F, cluster_columns = F)
dev.off()
Openfile("Evaluate.celltype.scores.random.pdf")
##############################################################################################################
#Extended data Figure 2d
require(GSVA)
Celltype.averages <- assays(Summary.DE.celltype)[["E"]]
APOE.paths.celltype.activity <- gsva(Celltype.averages, gset.idx.list = All.paths)
Celltype.assign <- colnames(APOE.paths.celltype.activity)[apply(APOE.paths.celltype.activity, 1, which.max)]
pdf("APOErelatedPathwayActivity.pdf", width = 2.5)
Heatmap(APOE.paths.celltype.activity[,c("Ex", "In", "Ast", "Mic", "Oli", "Opc")], show_row_names = F, name="pathway\nactivity", split = factor(Celltype.assign, levels=colnames(APOE.paths.celltype.activity)), cluster_columns = F, cluster_rows = T)
dev.off()
Openfile("APOErelatedPathwayActivity.pdf")
##############################################################################################################
#Extended data Figure 2e
set.seed(123)
APOE.paths.celltype.activity.random <- gsva(Celltype.averages, gset.idx.list = lapply(All.paths, function(i) sample(rownames(sce), length(i))))
pdf("APOErelatedPathwayActivity.random.pdf", width = 2.5)
Heatmap(APOE.paths.celltype.activity.random, show_row_names = F, name="pathway\nactivity")
dev.off()
Openfile("APOErelatedPathwayActivity.random.pdf")

pdf("APOErelatedPathwayActivity.boxplot.pdf", width = 3, height = 4)
boxplot(APOE.paths.celltype.activity, col=All.colors[colnames(APOE.paths.celltype.activity)], pch=20, las=2, ylab="pathway activity")
abline(h=0, lty=2)
boxplot(APOE.paths.celltype.activity.random, col=All.colors[colnames(APOE.paths.celltype.activity)], pch=20, las=2, ylab="pathway activity", main="random")
abline(h=0, lty=2)
dev.off()
Openfile("APOErelatedPathwayActivity.boxplot.pdf")
##############################################################################################################
