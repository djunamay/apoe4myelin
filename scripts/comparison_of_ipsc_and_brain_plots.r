########## analysis related to extended data figure 6 #############
###################################################################
print('|| comparison of ipsc celltypes and human brain - plotting... ||')

# required packages
library(ggpubr)
library(GSVA)
library(ComplexHeatmap)
library(ggplot2)
library(reshape2)

merged_scaled = read.csv('../data/supplementary_tables/ipsc_postmortem_merged_scaled_matrices_individual_level.csv', row.names = 'Row.names')
merged_scaled$X = NULL

# get colors and names
all.names = c()
for (i in 1:length(colnames(merged_scaled))){
    all.names = c(all.names, unlist(strsplit(colnames(merged_scaled)[i], "[.]"))[1])
}
model = unlist(lapply(all.names, function(x) strsplit(x,'_')[[1]][2]))
names = unlist(lapply(all.names, function(x) strsplit(x,'_')[[1]][1]))
names[names=='Neu'] = 'Ex'
col = readRDS('../data/single_cell_data/Cell_group_colors.rds')

####### do PC decomp
PCA <- prcomp(t(merged_scaled))
d = as.data.frame(PCA$x)

p <- ggplot(d,
  aes(x=PC1, y=PC2)) +
  geom_point(alpha = .5, aes(shape=model, color = names),size=5) + scale_colour_manual(values = (col[names]))

pdf('../plots/Extended_6/pca_ipsc.pdf', width = 5, height = 4)
p + theme_bw()  + theme(panel.background = element_rect(colour = "black", size=1), panel.grid.minor = element_blank(), panel.grid.major = element_blank())
dev.off()

write.csv(as.data.frame(summary(PCA)$importance[2,c('PC1','PC2')]), '../data/supplementary_tables/ipsc_post_mortem_pca_var_explained.csv')

####### distance in full gene space
distance = dist(t(merged_scaled), method = 'euclidean')

# get indices for human/ipsc and subset matrix
names = colnames(merged_scaled)
names = unlist(lapply(names, function(x) strsplit(x,'_')[[1]][2]))
h_index = unlist(lapply(names, function(x) startsWith(x, 'human')))
distance = as.matrix(distance)
new_dist = distance[!h_index,h_index]

# get cell order
i_cell_index = unlist(lapply(rownames(new_dist), function(x) strsplit(x,'_')[[1]][1]))
h_cell_index = unlist(lapply(colnames(new_dist), function(x) strsplit(x,'_')[[1]][1]))

# for each ipsc celltype, get the corresponding distances to the human celltypes
vals = list()
for(i in unique(i_cell_index)){
  index = which(i_cell_index==i)
  sele = new_dist[index,]
  for(x in unique(h_cell_index)){
    index = which(h_cell_index==x)
    curr = sele[,index]
    vals[[i]][[x]] = array(curr)
  }
}

# get distances of oligodendroglia to post-mortem cell types
d = do.call('rbind',vals$Oli)
f = melt(d)

# order celltypes by sum of ranks
f$rank = rank(f$value)
out = list()
for(i in unique(f$Var1)){
  out[[i]] = sum(f[f$Var1==i,'rank'])
}
out2 = unlist(out)
names(out2) = names(out)
order = out2[order(out2, decreasing=F)]

# invert values --> larger values indicate more similarity to oligodendroglia
f$value = max(f$value) - f$value
f$value = f$value/max(f$value) # show values as fraction of maximum distance from oligodendroglia
f$Var1 = factor(f$Var1, levels = names(order))

pdf('../plots/Extended_6/distplot.pdf', width = 3, height = 4)
ggplot(f, aes(x=Var1, y=value)) +
  geom_boxplot() + theme_classic() + stat_compare_means(method = 'wilcox.test', comparisons = list(c('Opc', 'Oli'), c('Opc', 'In'), c('Opc', 'Ex'), c('Opc', 'Mic'), c('Opc', 'Ast')))
dev.off()

####### cholesterol and myelination boxplot_w_stats
pathways = readRDS('../data/other_analyses_outputs/pathways.rds')

# look at the cholesterol biosynthesis signature
data = readRDS('../data/other_analyses_outputs/cholesterol_analysis.rds')
biosynth_genes = data[['union_cholest_biosynth']]$genes
myelination_genes = c('MYRF', 'MOG', 'PLP1', 'PLLP', 'MAG', 'OPALIN')

paths = list()
paths$cholest = biosynth_genes
paths$myelination = myelination_genes

# combine all the individual-level - celltype level gene expression profiles
gsva_out = gsva(as.matrix(merged_scaled), paths[c('cholest','myelination')], mx.diff=TRUE, verbose=TRUE, kcdf=c("Gaussian"), min.sz=5)
names = unlist(strsplit(colnames(gsva_out),'_')) %>% .[c(TRUE,FALSE)]

d = as.data.frame(t(gsva_out))
names[names == 'Neu'] = 'Ex'
d$celltype = names
d$model = model

order = c('Oli','Opc','Ast','Mic','Ex','In')

d$celltype = factor(d$celltype, levels = order)
d$grp = paste0(d$celltype, '_', d$model)
d$grp = factor(d$grp, levels = c('Oli_human','Oli_iPSC','Opc_human','Ast_human','Ast_iPSC','Mic_human','Mic_iPSC','Ex_human','Ex_iPSC','In_human','In_iPSC'))

col1 = c("#F39B7FFF" ,"#F39B7FFF" , "#8491B4FF", "#E64B35FF" , "#E64B35FF" , "#3C5488FF", "#3C5488FF", "#4DBBD5FF", "#4DBBD5FF","#00A087FF","#00A087FF")
names(col1) = levels(d$grp)

pdf('../plots/Extended_6/boxplot_cholest.pdf', width = 3,  height = 4)
ggplot(d, aes(x=grp, y=cholest)) +
  geom_boxplot() + theme_classic() + stat_compare_means(method = 'wilcox.test', comparisons = list(c('Oli_iPSC','Ast_iPSC'), c('Oli_iPSC', 'Mic_iPSC'), c('Oli_iPSC', 'Ex_iPSC')))
dev.off()

pdf('../plots/Extended_6/boxplot_myelin.pdf', width = 3,  height = 4)
ggplot(d, aes(x=grp, y=myelination)) +
  geom_boxplot() + theme_classic() + stat_compare_means(method = 'wilcox.test', comparisons = list(c('Oli_iPSC','Ast_iPSC'), c('Oli_iPSC', 'Mic_iPSC'), c('Oli_iPSC', 'Ex_iPSC')))
dev.off()

####### heatmaps
# get individual-level averages
names = c( 'Ast_iPSC','Mic_iPSC',  'Oli_iPSC', 'Neu_iPSC','Ast_human', 'Mic_human', 'Opc_human', 'Oli_human','In_human', 'Ex_human')
out = lapply(names, function(x) rowMeans(merged_scaled[,startsWith(colnames(merged_scaled), x)]))
names(out) = names
all_avs = t(do.call('rbind', out))
group = c('Ast','Mic','Oli','Ex','Ast','Mic', 'Opc', 'Oli', 'Ex', 'Ex')
c = c('grey','forestgreen')
names(c) = c('iPSC', 'human')
group2 = c('iPSC','iPSC','iPSC','iPSC','human','human','human','human','human','human')
column_ha = HeatmapAnnotation(group = group, group2= group2, col = list(group = col[group], group2 = c[group2]))
h1 = Heatmap(all_avs[paths$cholest,], border = T, rect_gp = gpar(col = 'black', lwd = 1), column_title = 'cholesterol biosynthesis',   column_title_gp = gpar(fontsize = 9),row_names_gp = gpar(fontsize = 11),column_names_gp = gpar(fontsize = 13), cluster_rows = T, bottom_annotation = column_ha)

pdf('../plots/Extended_6/cholest_genes_expression.pdf', width = 5, height = 5)
h1
dev.off()

h2 = Heatmap(all_avs[paths$myelination,], border = T, rect_gp = gpar(col = 'black', lwd = 1), column_title = 'myelination',   column_title_gp = gpar(fontsize = 9),row_names_gp = gpar(fontsize = 11),column_names_gp = gpar(fontsize = 13), cluster_rows = T, bottom_annotation = column_ha)

pdf('../plots/Extended_6/myelin_expression.pdf', width = 5, height = 3.8)
h2
dev.off()

print('done.')
