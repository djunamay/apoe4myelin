########## plots for Figure 2 ##########
########################################

source('../functions/plotting.r')
source('../functions/pathway_analyses.r')

##### required libraries ####
library('readxl')
library('ComplexHeatmap')
library('circlize')
library('tidyr')
library('ggplot2')
library('ggpubr')

print('loading data..')
data = readRDS('../data/other_analyses_outputs/cholesterol_analysis.rds')
cholest_genes_dist = data[['union_cholest_biosynth']]$density_plot_input
x = data$deg_level_analysis$degs

print('drawing the density plot')
# draw the density plot
pdf('../plots/Figure_2/cholest_biosynth_activity.pdf', width = 5 , height = 2.5)
ggplot(cholest_genes_dist, aes(x = biosynth)) + geom_density(aes(fill = factor(apoe_geno), alpha = .5))+ facet_wrap( ~ apoe_geno, nrow = 3) + geom_boxplot(width = .4) + geom_point(aes(x = biosynth,y=0, col = (AD))) + theme_classic() +
xlim(-1.5,1.5)
dev.off()

# save cor.test info to text file

temp = cor.test(cholest_genes_dist$biosynth, ifelse(cholest_genes_dist$apoe_geno == 33, 0, ifelse(cholest_genes_dist$apoe_geno == 34, 1, 2)))
write.csv(as.data.frame(do.call('rbind',(temp))), '../data/supplementary_tables/cholesterol_correlation_results.csv')

print('drawing the deg barplot')
# plotting the degs
pdf('../plots/Figure_2/cholesterol_degs.pdf', width = 3, height = 4)
barplot(x[order(x)], horiz = T, las = 1)
dev.off()

print('done.')
