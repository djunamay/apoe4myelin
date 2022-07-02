########## analysis related to extended data figure 6 #############
###################################################################
print('|| comparison of ipsc celltypes and human brain - get scaled matrices... ||')

# required packages
library(tidyr)
library(tibble)
library(reshape2)

# load all the data
human <- readRDS('../data/single_cell_data/individual_level_averages_per_celltype.rds')
ast <- as.data.frame(read.csv('../data/iPSC_data/FPKM_table_AST.txt', sep = '\t')) %>% column_to_rownames(., 'gene')
mic <- as.data.frame(read.csv('../data/iPSC_data/FPKM_table_MIC.txt', sep = '\t')) %>% column_to_rownames(., 'gene')
oligodendroglia <- as.data.frame(read.csv('../data/iPSC_data/FPKM_table_OPC.txt', sep = '\t')) %>% column_to_rownames(., 'gene')
neuron <- as.data.frame(read.csv('../data/iPSC_data/FPKM_table_NEU.txt', sep = '\t'))

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

write.csv(merged_scaled, '../data/supplementary_tables/ipsc_postmortem_merged_scaled_matrices_individual_level.csv')
print('done.')
