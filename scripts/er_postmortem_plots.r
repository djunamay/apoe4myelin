########## ER stress related to extended data figure 8 #############
##########################################################################
print('|| plotting ER stress results... ||')

source('../functions/plotting.r')

##### required libraries ####
library('readxl')
library('ComplexHeatmap')
library('circlize')
library('tidyr')
library('ggplot2')
library('ggpubr')

data = readRDS('../data/other_analyses_outputs/er_stress_results.rds')

# plot the ATF6 pathways
df = data$ER_stress$ATF6_activity
pdf('../plots/Extended_7/unfolded_protein_response.pdf', width = 2, height = 3)
p <- ggplot(df, aes(x=APOE, y=value)) +
  geom_boxplot()
p + geom_jitter(shape=16, position=position_jitter(0.2)) + stat_compare_means(method = "wilcox.test")
dev.off()

x = data$ER_stress$ATF6_pathway_degs
pdf('../plots/Extended_7/unfolded_genes.pdf', width = 3.5, height = 4)
barplot(x[order(x)], horiz = T, las = 1)
dev.off()

print('done')
