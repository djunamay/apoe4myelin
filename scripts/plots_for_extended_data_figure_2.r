########## plots for extended data figure 2 ##########
######################################################

library(SingleCellExperiment)
library(ComplexHeatmap)
library(GSVA)

source('../functions/qc_and_annotation_aux_functions.r')

#Evaluate celltype scores
Celltype.averages <- readRDS('../data/single_cell_data/individual_level_averages_per_celltype.rds')
temp = lapply(names(Celltype.averages), function(x) rowMeans(Celltype.averages[[x]]))
names(temp) = names(Celltype.averages)
Celltype.averages = do.call('cbind', temp)
RefCellTypeMarkers = readRDS('../data/single_cell_data/RefCellTypeMarkers.adultBrain.rds')
sce = readRDS('../data/single_cell_data/single_cell_experiment_object_qced.rds')

pdf("../plots/Extended_2/Evaluate.celltype.scores.pdf", height = 2, width = 3)
Celltype.activity <- t(gsva(Celltype.averages, gset.idx.list = RefCellTypeMarkers))
Celltype.activity <- Celltype.activity[colnames(Celltype.activity),]
Heatmap(Celltype.activity, cluster_rows = F, cluster_columns = F, rect_gp = gpar(col = "black"), name="geneset\nactivity score")
Map <- Table.to.matrix(table(colnames(Celltype.activity)[apply(Celltype.activity, 1, which.max)], rownames(Celltype.activity)))
Heatmap(Map, col=c("white","black"), rect_gp = gpar(col = "black"), name="assignment", cluster_rows = F, cluster_columns = F)
dev.off()

#Evaluate celltype scores random
#   subscript out of bounds
pdf("../plots/Extended_2/Evaluate.celltype.scores.random.pdf", height = 2, width = 3)
Celltype.activity.random <- t(gsva(Celltype.averages, gset.idx.list = lapply(RefCellTypeMarkers, function(i) sample(rownames(sce), length(i)))))
Celltype.activity.random <- Celltype.activity.random[colnames(Celltype.activity),]
Heatmap(Celltype.activity.random[colnames(Celltype.activity), rownames(Celltype.activity)], rect_gp = gpar(col = "black"), cluster_rows = F, cluster_columns = F, name="geneset\nactivity score")
Maprandom <- Table.to.matrix(table(colnames(Celltype.activity.random)[apply(Celltype.activity.random, 1, which.max)], rownames(Celltype.activity.random)))
Heatmap(Maprandom, col=c("white","black"), rect_gp = gpar(col = "black"), name="assignment", cluster_rows = F, cluster_columns = F)
dev.off()

#APOE_related_pathway activity
require(GSVA)
pathways = readRDS('../data/other_analyses_outputs/pathways.rds')
All.paths = pathways$pathways$apoe_gsets_all

Celltype.averages <- assays(Summary.DE.celltype)[["E"]]
APOE.paths.celltype.activity <- gsva(Celltype.averages, gset.idx.list = All.paths)
Celltype.assign <- colnames(APOE.paths.celltype.activity[,c("Ex", "In", "Ast", "Mic", "Oli", "Opc")])[apply(APOE.paths.celltype.activity[,c("Ex", "In", "Ast", "Mic", "Oli", "Opc")], 1, which.max)]
pdf("../plots/Extended_2/APOErelatedPathwayActivity.pdf", width = 2.5)
Heatmap(APOE.paths.celltype.activity[,c("Ex", "In", "Ast", "Mic", "Oli", "Opc")], show_row_names = F, name="pathway\nactivity", split = factor(Celltype.assign, levels=colnames(APOE.paths.celltype.activity[,c("Ex", "In", "Ast", "Mic", "Oli", "Opc")])), cluster_columns = F, cluster_rows = T)
dev.off()

#APOE-related pathway activity random
APOE.paths.celltype.activity.random <- gsva(Celltype.averages, gset.idx.list = lapply(All.paths, function(i) sample(rownames(sce), length(i))))
pdf("../plots/Extended_2/APOErelatedPathwayActivity.random.pdf", width = 2.5)
Heatmap(APOE.paths.celltype.activity.random[,c("Ex", "In", "Ast", "Mic", "Oli", "Opc")], show_row_names = F, name="pathway\nactivity")
dev.off()

# apoe-related pathway activity boxplots
All.colors = readRDS('../data/single_cell_data/Cell_group_colors.rds')

pdf("../plots/Extended_2/APOErelatedPathwayActivity.boxplot.pdf", width = 3, height = 4)
boxplot(APOE.paths.celltype.activity[,c("Ex", "In", "Ast", "Mic", "Oli", "Opc")], col=All.colors[colnames(APOE.paths.celltype.activity)], pch=20, las=2, ylab="pathway activity")
abline(h=0, lty=2)
boxplot(APOE.paths.celltype.activity.random[,c("Ex", "In", "Ast", "Mic", "Oli", "Opc")], col=All.colors[colnames(APOE.paths.celltype.activity)], pch=20, las=2, ylab="pathway activity", main="random")
abline(h=0, lty=2)
dev.off()
##############################################################################################################
