########## plots for extended data figure 2 ##########
######################################################

#Evaluate celltype scores
Summary.DE.celltype <- readRDS("./share_data/Summary.DE.celltype.rds")
Celltype.averages <- assays(Summary.DE.celltype)[["E"]]

pdf("../plots/Extended_2/Evaluate.celltype.scores.pdf", height = 2, width = 3)
Celltype.activity <- gsva(Celltype.averages, gset.idx.list = RefCellTypeMarkers)
Celltype.activity <- Celltype.activity[colnames(Celltype.activity),]
Heatmap(Celltype.activity, cluster_rows = F, cluster_columns = F, rect_gp = gpar(col = "black"), name="geneset\nactivity score")
Map <- Table.to.matrix(table(colnames(Celltype.activity)[apply(Celltype.activity, 1, which.max)], rownames(Celltype.activity)))
Heatmap(Map, col=c("white","black"), rect_gp = gpar(col = "black"), name="assignment", cluster_rows = F, cluster_columns = F)
dev.off()

#Evaluate celltype scores random
pdf("../plots/Extended_2/Evaluate.celltype.scores.random.pdf", height = 2, width = 3)
Celltype.activity.random <- gsva(Celltype.averages, gset.idx.list = lapply(RefCellTypeMarkers, function(i) sample(rownames(sce), length(i))))
Celltype.activity.random <- Celltype.activity.random[colnames(Celltype.activity),]
Heatmap(Celltype.activity.random[colnames(Celltype.activity), rownames(Celltype.activity)], rect_gp = gpar(col = "black"), cluster_rows = F, cluster_columns = F, name="geneset\nactivity score")
Maprandom <- Table.to.matrix(table(colnames(Celltype.activity.random)[apply(Celltype.activity.random, 1, which.max)], rownames(Celltype.activity.random)))
Heatmap(Maprandom, col=c("white","black"), rect_gp = gpar(col = "black"), name="assignment", cluster_rows = F, cluster_columns = F)
dev.off()

#APOE_related_pathway activity
require(GSVA)
Celltype.averages <- assays(Summary.DE.celltype)[["E"]]
APOE.paths.celltype.activity <- gsva(Celltype.averages, gset.idx.list = All.paths)
Celltype.assign <- colnames(APOE.paths.celltype.activity)[apply(APOE.paths.celltype.activity, 1, which.max)]
pdf("../plots/Extended_2/APOErelatedPathwayActivity.pdf", width = 2.5)
Heatmap(APOE.paths.celltype.activity[,c("Ex", "In", "Ast", "Mic", "Oli", "Opc")], show_row_names = F, name="pathway\nactivity", split = factor(Celltype.assign, levels=colnames(APOE.paths.celltype.activity)), cluster_columns = F, cluster_rows = T)
dev.off()

#APOE-related pathway activity random
set.seed(123)
APOE.paths.celltype.activity.random <- gsva(Celltype.averages, gset.idx.list = lapply(All.paths, function(i) sample(rownames(sce), length(i))))
pdf("../plots/Extended_2/APOErelatedPathwayActivity.random.pdf", width = 2.5)
Heatmap(APOE.paths.celltype.activity.random, show_row_names = F, name="pathway\nactivity")
dev.off()
Openfile("APOErelatedPathwayActivity.random.pdf")

# apoe-related pathway activity boxplots
pdf("../plots/Extended_2/APOErelatedPathwayActivity.boxplot.pdf", width = 3, height = 4)
boxplot(APOE.paths.celltype.activity, col=All.colors[colnames(APOE.paths.celltype.activity)], pch=20, las=2, ylab="pathway activity")
abline(h=0, lty=2)
boxplot(APOE.paths.celltype.activity.random, col=All.colors[colnames(APOE.paths.celltype.activity)], pch=20, las=2, ylab="pathway activity", main="random")
abline(h=0, lty=2)
dev.off()
##############################################################################################################
