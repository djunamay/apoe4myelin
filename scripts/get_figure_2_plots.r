########## script 5 in APOE4 impairs myelination via cholesterol dysregulation in oligodendrocytes ##########
#############################################################################################################

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
data = readRDS('../data/cholesterol_analysis.rds')
cholest_genes_dist = data[['union_cholest_biosynth']]$density_plot_input
x = data$deg_level_analysis$degs

print('drawing the density plot')
# draw the density plot
pdf('../plots/cholest_biosynth_activity.pdf', width = 5 , height = 2.5)
ggplot(cholest_genes_dist, aes(x = biosynth)) + geom_density(aes(fill = factor(apoe_geno), alpha = .5))+ facet_wrap( ~ apoe_geno, nrow = 3) + geom_boxplot(width = .4) + geom_point(aes(x = biosynth,y=0, col = (AD))) + theme_classic() +
xlim(-1.5,1.5)
dev.off()

# save cor.test info to text file
cor.test(cholest_genes_dist$biosynth, ifelse(cholest_genes_dist$apoe_geno == 33, 0, ifelse(cholest_genes_dist$apoe_geno == 34, 1, 2)))

print('drawing the deg barplot')
# plotting the degs
pdf('../plots/cholesterol_degs.pdf', width = 3, height = 4)
barplot(x[order(x)], horiz = T, las = 1)
dev.off()

# plot the ATF6 pathways
df = data$ER_stress$ATF6_activity
pdf('../plots/unfolded_protein_response.pdf', width = 2, height = 3)
p <- ggplot(df, aes(x=APOE, y=value)) +
  geom_boxplot()
p + geom_jitter(shape=16, position=position_jitter(0.2)) + stat_compare_means(method = "wilcox.test")
dev.off()
x = data$ER_stress$ATF6_pathway_degs
pdf('../plots/unfolded_genes.pdf', width = 3.5, height = 4)
barplot(x[order(x)], horiz = T, las = 1)
dev.off()

print('done')
