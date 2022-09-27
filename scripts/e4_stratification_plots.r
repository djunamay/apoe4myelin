########## plots of effects stratification for extended data figure 3 #############
##########################################################################

print('|| plotting stratification results. ||')

library(ggplot2)

all_data = readRDS('../data/other_analyses_outputs/stratified_anaylsis.rds')

# show the APOE4 effects on nonAD background as barplot
df = all_data[['res']][['lipid_associated']][['APOE34_effect_nia34']]
df = df[df$P.Value<0.05,]
x = (df$logFC)
names(x) = rownames(df)
pdf('../plots/Extended_2/APOE34_lipid_effects_no_path.pdf', width = 3, height = 6)
barplot(x[order(x)], horiz = T, las = 1)
dev.off()

# show boxplots for the cholesterol pathway of interest
out_lipid_terms = all_data[['res']][['lipid_associated']][['gsva_out']]
df = as.data.frame(out_lipid_terms$Oli[,'cholesterol biosynthesis III (via desmosterol) Homo sapiens PWY66-4'])
summary = read.csv('../data/single_cell_data/metadata_by_individual.csv')
rownames(summary) = summary[,'projid']

colnames(df) = c('pathway')
df$apoe_genotype = summary[rownames(df), 'apoe_genotype']
df$APOE = ifelse(df$apoe_genotype%in%c(33), 'E3', 'E4')
df$AD = (summary[rownames(df), 'AD'])
df$AD = ifelse(df$AD%in%c('yes'),'AD','non AD')
df$AD = factor(df$AD, levels = c('non AD','AD'))
df$both = paste0(df$APOE, '_', df$AD)
df$both = factor(df$both, levels = c('E3_non AD', 'E3_AD', 'E4_non AD', 'E4_AD'))

x = list()
x[['p_APOE34_effect_AD']] <- ggplot(df[df$AD=='AD' & df$apoe_genotype!=44,], aes(x=factor(APOE), y=pathway)) +
  geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_classic()

x[['p_APOE34_effect_noAD']] <- ggplot(df[df$AD=='non AD'& df$apoe_genotype!=44,], aes(x=factor(APOE), y=pathway)) +
  geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_classic()

x[['p_AD_effect_unstratified']] <- ggplot(df, aes(x=factor(AD), y=pathway)) +
  geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_classic()

x[['p_AD_effect_APOE34']] <- ggplot(df[df$apoe_genotype==34,], aes(x=factor(AD), y=pathway)) +
  geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_classic()#+ stat_compare_means('wilcox.test')

x[['p_AD_effect_APOE33']] <- ggplot(df[df$APOE=='E3',], aes(x=factor(AD), y=pathway)) +
  geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_classic()#+ stat_compare_means('wilcox.test')

x[['p_AD_and_APOE34_effects']] <- ggplot(df[df$apoe_genotype!=44,], aes(x=factor(both), y=pathway)) +
  geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_classic()#+ stat_compare_means('wilcox.test')

for(i in names(x)){
    pdf(paste0('../plots/Extended_2/',i,'.pdf'), width = 2, height = 2)
    print(x[[i]])
    dev.off()
}

print('done.')
