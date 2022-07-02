########## iPSC analysis related to extended data figure 8 #############
##########################################################################
print('|| plotting iPSC deg results of interest... ||')

# required packages
library(GSVA)
library(limma)
library(tidyr)
library(reshape2)
library(ggplot2)

out_APOE4 = read.table('../data/iPSC_data/OPC_DEG_statistics.txt',  header = T)
pathways = readRDS('../data/other_analyses_outputs/pathways.rds')

# look at the cholesterol biosynthesis signature
data = readRDS('../data/other_analyses_outputs/cholesterol_analysis.rds')
biosynth_genes = data[['union_cholest_biosynth']]$genes

out = list()
out[['cholest_bio']] = biosynth_genes
out[['ATF6']] = c('HSP90B1', 'CALR', 'HSPA5', 'MBTPS1', 'ATF6')

curr = out_APOE4
curr = curr[!duplicated(curr$gene_id),]
rownames(curr) = curr$gene_id

l = list()
for(i in names(out)){
    df = na.omit(curr[out[[i]],])
    df$score = df$log2.fold_change.
    df$grp = i
    l[[i]] = df[df$q_value<0.05,]
}

for(f in names(l)){
    x = l[[f]]$score
    names(x) = rownames(l[[f]])
    pdf(paste0('../plots/Extended_8/', f,'.pdf'), width = 3, height = 5)
    print(barplot(x[order(x)], las = 1, horiz = T))
    dev.off()
}

# show other key genes of interest
norm_counts = read.csv('../data/iPSC_data/FPKM_table_OPC.txt', sep = '\t')
rownames(norm_counts) = norm_counts$gene
norm_counts$gene = NULL
df = na.omit(norm_counts[c('SOAT1', 'SOAT2', 'CYP46A1'),])
df$gene = rownames(df)
x = melt(df)
x$variable = ifelse(startsWith(as.character(x$variable), 'E3'), 'APOE3', 'APOE4')

pdf('../plots/Extended_8/SOAT1_CYP_boxplots_ipsc.pdf', width = 3, height = 4)
ggplot(x, aes(x=variable, y=value, col = variable)) +
      geom_boxplot(width = .5) + geom_jitter() + facet_wrap(. ~ gene,  scales="free_y", nrow = 2) + theme_classic()
dev.off()

write.csv(out_APOE4[out_APOE4$gene%in%c('SOAT1', 'CYP46A1'),], '../data/supplementary_tables/SOAT1_CYP_stats_ipsc.csv')

print('done.')
