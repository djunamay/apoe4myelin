##### PCA ######

library(tidyr)
library(tibble)
library(ggplot2)

# load all the data
human <- readRDS('../../submission_code_09012021/data/Averages.by.celltype.by.individual.rds')
ast <- as.data.frame(read.csv('../../submission_code_09012021/data/FPKM_table_AST.txt', sep = '\t')) %>% column_to_rownames(., 'gene')
mic <- as.data.frame(read.csv('../../submission_code_09012021/data/FPKM_table_MIC.txt', sep = '\t')) %>% column_to_rownames(., 'gene')
oligodendroglia <- as.data.frame(read.csv('../../submission_code_09012021/data/FPKM_table_OPC.txt', sep = '\t')) %>% column_to_rownames(., 'gene')
neuron <- as.data.frame(read.csv('../../submission_code_09012021/data/DEG_CTRL_E3_NEU_CTRL_E4_NEU_FPKM.txt', sep = '\t'))

# get mean expression profiles for human data
human = human[c('Ast', 'Mic', 'Opc', 'Oli', 'Ex', 'In')]
gene_names = Reduce(intersect, list(rownames(human$Ast), rownames(human$Mic), rownames(human$Opc), rownames(human$Oli)))
human_i = lapply(names(human), function(i) human[[i]][gene_names,])
names(human_i) = names(human)
d_i_h = as.data.frame(do.call(cbind, human_i))
h_names = c()
for (i in names(human_i)){
    h_names = c(h_names, rep(paste0(i,'_','human'), 32))
}
colnames(d_i_h) = h_names
h_scaled = na.omit(t(scale(t(d_i_h))))

# compile the ipsc data
av_i = list()
av_i[['Ast.ipsc' ]] = ast
av_i[['Mic.ipsc' ]] = mic
av_i[['Oli.ipsc' ]] = oligodendroglia
av_i[['Neu.ipsc' ]] = neuron
gene_names = Reduce(intersect, list(rownames(av_i$Ast.ipsc), rownames(av_i$Mic.ipsc), rownames(av_i$Oli.ipsc),  rownames(av_i$Neu.ipsc)))
av_ii = lapply(names(av_i), function(i) av_i[[i]][gene_names,])
names(av_ii) = names(av_i)
d_i_i = as.data.frame(do.call(cbind, av_ii))
new.names = c()
names = strsplit(colnames(d_i_i),'[.]')
for (i in 1:length(names)){
    new.names = c(new.names,paste0(names[[i]][1], '_','iPSC'))
}
colnames(d_i_i) = new.names
i_scaled = na.omit(t(scale(t(d_i_i))))

# merge the scaled matrices
merged_scaled = merge(i_scaled,h_scaled, by = 0, how = 'inner')
merged_scaled$Row.names = NULL

# get colors and names
all.names = c()
for (i in 1:length(colnames(merged_scaled))){
    all.names = c(all.names, unlist(strsplit(colnames(merged_scaled)[i], "[.]"))[1])
}

model = unlist(lapply(all.names, function(x) strsplit(x,'_')[[1]][2]))
names = unlist(lapply(all.names, function(x) strsplit(x,'_')[[1]][1]))
names[names=='Neu'] = 'Ex'
col = readRDS('../../submission_code_09012021/data/Cell_group_colors.rds')

# do PC decomp
PCA <- prcomp(t(merged_scaled))
d = as.data.frame(PCA$x)

p <- ggplot(d,
  aes(x=PC1, y=PC2)) +
  geom_point(size=3, alpha = .5, aes(shape=model, color = names),size=5) + scale_colour_manual(values = (col[names]))

pdf('../plots/pca_ipsc.pdf', width = 5, height = 4)
p + theme_bw()  + theme(panel.background = element_rect(colour = "black", size=1), panel.grid.minor = element_blank(), panel.grid.major = element_blank())
dev.off()

p <- ggplot(d, aes(x=PC1, y=PC2, shape = model)) +
  geom_point(size=3, alpha = 1, aes(color = names),size=10) +
  geom_point(colour = "white",alpha = .5, size = 1.5) + scale_colour_manual(values = (col[names]))+ theme_bw()  + theme(panel.background = element_rect(colour = "black", size=1), panel.grid.minor = element_blank(), panel.grid.major = element_blank())

print('summary PCs:')
print(summary(PCA)$importance[2,c('PC1','PC2')])

###### Distances in Full gene space ########
# compute distances in gene space
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

library(reshape2)
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

# normalize so max distance is 1
d = do.call('rbind',vals$Oli)
f = melt(d)

# order by sum of ranks (wilcox statistic)
f$rank = rank(f$value)
out = list()
for(i in unique(f$Var1)){
  out[[i]] = sum(f[f$Var1==i,'rank'])
}
out2 = unlist(out)
names(out2) = names(out)

order = out2[order(out2, decreasing=F)]

f$value = max(f$value) - f$value
f$value = f$value/max(f$value)
f$Var1 = factor(f$Var1, levels = names(order))

# sorted by means
pdf('../plots/distplot.pdf', width = 3, height = 4)
boxplot((value) ~ Var1,f, col = col[levels(f$X1)], ylab = 'd(Oli ipsc, X) in gene space', xlab = 'x', las = 2)
dev.off()

# get corresponding wilcox p-values
wilcox.test(vals$Oli$Opc, vals$Oli$Oli, conf.int = T)

wilcox.test(vals$Oli$Opc, vals$Oli$In, conf.int = T)

wilcox.test(vals$Oli$Opc, vals$Oli$Ast, conf.int = T)

wilcox.test(vals$Oli$Opc, vals$Oli$Ex, conf.int = T)

wilcox.test(vals$Oli$Opc, vals$Oli$Mic, conf.int = T)


######## cholesterol and myelination boxplot_w_stats
pathways = readRDS('../data_outputs/pathways.rds')

# look at the cholesterol biosynthesis signature
f1_data = readRDS('../data_outputs/pathway_scores.rds')
oli.lipid.fits = f1_data$res$lipid_associated$fits$Oli
lipid.paths.oli = pathways$pathways$low_removed_lipid_associated$Oli
paths = rownames(oli.lipid.fits[oli.lipid.fits$P.Value<0.05,])
biosynth_genes = unique(unname(unlist(lipid.paths.oli[paths])))
myelination_genes = c('MYRF', 'MOG', 'PLP1', 'PLLP', 'MAG', 'OPALIN')

paths = list()
paths$cholest = biosynth_genes
paths$myelination = myelination_genes

# combine all the individual-level - celltype level gene expression profiles
library(GSVA)
merged_scaled = merge(i_scaled,h_scaled, by = 0, how = 'inner')
rownames(merged_scaled) = merged_scaled$Row.names
merged_scaled$Row.names = NULL

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

pdf('../plots/boxplot_cholest.pdf', width = 3,  height = 4)
boxplot(cholest ~ grp, d, col = col1[levels(d$grp)], ylab = 'gene set activity score', las = 2, main = 'cholesterol biosynthesis')
dev.off()

pdf('../plots/boxplot_myelin.pdf', width = 3,  height = 4)
boxplot(myelination ~ grp, d, col = col1[levels(d$grp)], ylab = 'gene set activity score', las = 2, main = 'myelination')
dev.off()

wilcox.test(d[d$grp=='Oli_iPSC','cholest'], d[d$grp=='Opc_human','cholest'])
wilcox.test(d[d$grp=='Oli_iPSC','cholest'], d[d$grp=='Oli_human','cholest'])

wilcox.test(d[d$grp=='Oli_iPSC','myelination'], d[d$grp=='Opc_human','myelination'])
wilcox.test(d[d$grp=='Oli_iPSC','myelination'], d[d$grp=='Oli_human','myelination'])

############ heatmaps
library(ComplexHeatmap)
#show cholest biosynthesis and myelination activity heatmaps on average values}
# get mean expression profiles for both systems data
human = human[c('Ast', 'Mic', 'Opc', 'Oli', 'Ex', 'In')]
gene_names = Reduce(intersect, list(rownames(human$Ast), rownames(human$Mic), rownames(human$Opc), rownames(human$Oli)))
human_i = lapply(names(human), function(i) human[[i]][gene_names,])
names(human_i) = names(human)

av_h = lapply(names(human_i), function(x) rowMeans(human_i[[x]]))
names(av_h) = names(human_i)
av_h = as.data.frame(do.call(cbind, av_h))
colnames(av_h) = c('Ast.human', 'Mic.human', 'Opc.human', 'Oli.human', 'Ex.human', 'In.human')
d_human = na.omit(t(scale(t(av_h))))

av_i = list()
av_i[['Ast.ipsc' ]] = ast
av_i[['Mic.ipsc' ]] = mic
av_i[['Oli.ipsc' ]] = oligodendroglia
av_i[['Neu.ipsc' ]] = neuron

gene_names = Reduce(intersect, list(rownames(av_i$Ast.ipsc), rownames(av_i$Mic.ipsc), rownames(av_i$Oli.ipsc),  rownames(av_i$Neu.ipsc)))
av_ii = lapply(names(av_i), function(i) av_i[[i]][gene_names,])
names(av_ii) = names(av_i)

d_i = lapply(names(av_ii), function(x) rowMeans(av_ii[[x]]))
names(d_i) = names(av_ii)
d_i = as.data.frame(do.call(cbind, d_i))
d_ipsc = na.omit(t(scale(t(d_i))))

df = merge(d_ipsc, d_human, by = 0, how = 'inner') %>% column_to_rownames(.,'Row.names')

group = c('Ast','Mic','Oli','Ex','Ast','Mic', 'Opc', 'Oli', 'Ex', 'Ex')
c = c('grey','forestgreen')
names(c) = c('iPSC', 'human')
group2 = c('iPSC','iPSC','iPSC','iPSC','human','human','human','human','human','human')
column_ha = HeatmapAnnotation(group = group, group2= group2, col = list(group = col[group], group2 = c[group2]))
h1 = Heatmap(df[paths$cholest,], border = T, rect_gp = gpar(col = 'black', lwd = 1), column_title = 'cholesterol biosynthesis',   column_title_gp = gpar(fontsize = 9),row_names_gp = gpar(fontsize = 11),column_names_gp = gpar(fontsize = 13), cluster_rows = T, bottom_annotation = column_ha)

pdf('../plots/cholest_genes_expression.pdf', width = 5, height = 5)
h1
dev.off()

h2 = Heatmap(df[paths$myelination,], border = T, rect_gp = gpar(col = 'black', lwd = 1), column_title = 'myelination',   column_title_gp = gpar(fontsize = 9),row_names_gp = gpar(fontsize = 11),column_names_gp = gpar(fontsize = 13), cluster_rows = T, bottom_annotation = column_ha)

pdf('../plots/myelin_expression.pdf', width = 5, height = 3.8)
h2
dev.off()
