########## lipidomic analysis for extended data figure 4 #############
#######################################################################
print('|| pfc_lipidomic_analysisis - plotting... ||')

source('../functions/plotting.r')

# required packages
library(readxl)
library(ggplot2)
library(ggpubr)

data_subset = read.csv('../data/supplementary_tables/pfc_lipidomics_data_metadata.csv')
df = read.csv('../data/supplementary_tables/pfc_lipidomics_data.csv')

# make metadata distribution plots
data_subset$APOE = ifelse(data_subset$apoe_genotype%in%c('33','23','22'), 'APOE4-noncarrier', 'APOE4-carrier')
data_subset$grp = paste0(data_subset$AD, '_', data_subset$APOE)
data_subset$grp = factor(data_subset$grp, levels = c('noAD_APOE4-noncarrier','AD_APOE4-noncarrier','noAD_APOE4-carrier','AD_APOE4-carrier'))
cols = c('grey', 'pink', 'black', 'red')
names(cols) =c('noAD_APOE4-noncarrier','AD_APOE4-noncarrier','noAD_APOE4-carrier','AD_APOE4-carrier')

pdf('../plots/Extended_4/pmi_pfc_lipid_cohort.pdf', width = )
boxplot_w_stats(df = data_subset, x = 'grp', y = 'pmi', palette = cols, comparisons = list(c('noAD_APOE4-noncarrier', 'noAD_APOE4-carrier'),c('AD_APOE4-noncarrier', 'AD_APOE4-carrier')), xlab = '', ylab = 'pmi', width = .5, alpha = .5)
dev.off()

pdf('../plots/Extended_4/nft_pfc_lipid_cohort.pdf', width = )
boxplot_w_stats(df = data_subset, x = 'grp', y = 'nft', palette = cols, comparisons = list(c('noAD_APOE4-noncarrier', 'noAD_APOE4-carrier'),c('AD_APOE4-noncarrier', 'AD_APOE4-carrier')), xlab = '', ylab = 'nft', width = .5, alpha = .5)
dev.off()

pdf('../plots/Extended_4/amyloid_pfc_lipid_cohort.pdf', width = )
boxplot_w_stats(df = data_subset, x = 'grp', y = 'amyloid', palette = cols, comparisons = list(c('noAD_APOE4-noncarrier', 'noAD_APOE4-carrier'),c('AD_APOE4-noncarrier', 'AD_APOE4-carrier')), xlab = '', ylab = 'amyloid', width = .5, alpha = .5)
dev.off()

pdf('../plots/Extended_4/amyloid_age_death_lipid_cohort.pdf', width = )
boxplot_w_stats(df = data_subset, x = 'grp', y = 'age_death', palette = cols, comparisons = list(c('noAD_APOE4-noncarrier', 'noAD_APOE4-carrier'),c('AD_APOE4-noncarrier', 'AD_APOE4-carrier')), xlab = '', ylab = 'age_death', width = .5, alpha = .5)
dev.off()

pdf('../plots/Extended_4/amyloid_msex_lipid_cohort.pdf')
get_barplot(data_subset, 'grp', 'msex')
dev.off()

# show the plots for lipid species of interest (pval < 0.1)
print(df[df$p.value<0.05,])

data_subset$APOE = factor(data_subset$APOE, levels = c('APOE4-noncarrier','APOE4-carrier'))
data_subset$AD = factor(data_subset$AD, levels = c('noAD', 'AD'))
pdf('../plots/Extended_4/ChE(18:1)_1.pdf', width = 10, height = 5)
ggplot(data_subset, aes(x=(APOE), y=(data_subset[,'ChE.18.1._1']), col = APOE)) +
stat_compare_means(comparisons = list(c('APOE4-noncarrier', 'APOE4-carrier')),method = "wilcox.test") + geom_violin() +  geom_boxplot(width = .1) + geom_jitter(size = .5) + theme_classic() + facet_grid(~ AD)
dev.off()

pdf('../plots/ChE.18.1._1_boxplot.pdf', width = 4, height = 7)
ggplot(data_subset, aes(x=(APOE), y=data_subset[,'ChE.18.1._1'], col = APOE)) +
stat_compare_means(comparisons = list(c('APOE4-noncarrier', 'APOE4-carrier')),method = "wilcox.test") +  geom_boxplot(width = 1) + theme_classic() + facet_grid(~ AD)
dev.off()

print(table(data_subset$APOE, data_subset$AD))

print(table(data_subset$apoe_genotype))

print(table(data_subset$niareagansc))
print('done.')
