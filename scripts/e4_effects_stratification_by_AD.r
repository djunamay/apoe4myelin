########## effects stratification for extended data figure 3 #############
##########################################################################
# load functions
source('../functions/pathway_analyses.r')
print('loading libraries..')

# required packages
library("readxl")
library('SingleCellExperiment')
library('GSVA')
library('limma')
library('parallel')
library('tidyr')
library(stringr)
library(GSA)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggpubr)

# load the data
print('loading data')
data = readRDS('../data/other_analyses_outputs/pathways.rds')
low_removed_bp = data$pathways$low_removed_bp
apoe_gsets_low_removed = data$pathways$apoe_gsets_low_removed
apoe_gsets = data[['pathways']][['apoe_gsets_all']]
all_paths = data$pathways$all

# output data
all_data = list()

# load the data
print('loading the data..')
expressed = readRDS('../data/single_cell_data/expressed_genes_per_celltype.rds')
av_expression = readRDS('../data/single_cell_data/individual_level_averages_per_celltype.rds')
summary = read.csv('../data/single_cell_data/metadata_by_individual.csv')
summary$APOE4 = ifelse(summary$apoe_genotype == '33','e3','e4')
rownames(summary) = summary[,'projid...2']
low_removed_lipid_paths = data$pathways$low_removed_lipid_associated

# run GSVA on lipid pathways pathways
print('running GSVA...')
order = c('Oli')
out_lipid_terms = lapply(order, function(x) t(gsva(as.matrix(av_expression[[x]]), low_removed_lipid_paths[[x]], mx.diff=TRUE, verbose=F, kcdf=c("Gaussian"), min.sz=5, max.sz = 150, parallel.sz=50)))
names(out_lipid_terms) = order
all_data[['res']][['lipid_associated']][['gsva_out']] = out_lipid_terms

# stratify by AD status
summary1 = summary[summary$niareagansc%in%c(3,4)& summary$apoe_genotype!=44,]
out_lipid_terms2 = lapply(names(out_lipid_terms), function(x) out_lipid_terms[[x]][rownames(out_lipid_terms[[x]])%in%as.character(rownames(summary1)),])
names(out_lipid_terms2) = names(out_lipid_terms)
fits = get_fits(out_lipid_terms2, summary1)
all_data[['res']][['lipid_associated']][['APOE34_effect_nia34']] = fits$Oli

summary2 = summary[summary$niareagansc%in%c(1,2) & summary$apoe_genotype!=44,]
out_lipid_terms2 = lapply(names(out_lipid_terms), function(x) out_lipid_terms[[x]][rownames(out_lipid_terms[[x]])%in%as.character(rownames(summary2)),])
names(out_lipid_terms2) = names(out_lipid_terms)
fits = get_fits(out_lipid_terms2, summary2)
all_data[['res']][['lipid_associated']][['APOE34_effect_nia12']] = fits$Oli

# show the APOE4 effects on nonAD background as barplot
df = all_data[['res']][['lipid_associated']][['APOE34_effect_nia34']]
df = df[df$P.Value<0.05,]
x = (df$logFC)
names(x) = rownames(df)
pdf('../plots/Extended_3/APOE34_lipid_effects_no_path.pdf', width = 3, height = 6)
barplot(x[order(x)], horiz = T, las = 1)
dev.off()

# show boxplots for the cholesterol pathway of interest
print('getting bar plots')
df = as.data.frame(out_lipid_terms$Oli[,'cholesterol biosynthesis III (via desmosterol) Homo sapiens PWY66-4'])
colnames(df) = c('pathway')
df$apoe_genotype = summary[rownames(df), 'apoe_genotype']
df$APOE = ifelse(df$apoe_genotype%in%c(33), 'E3', 'E4')
df$AD = (summary[rownames(df), 'niareagansc'])
df$AD = ifelse(df$AD%in%c(1,2),'AD','non AD')
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
    pdf(paste0('../plots/Extended_3/',i,'.pdf'), width = 2, height = 2)
    print(x[[i]])
    dev.off()
}

# get linear model fits for AD (AD effects on the background of APOE4 or E3)
print('getting linear model fits...')
summary$APOE4 = ifelse(summary$apoe_genotype == '33','e3','e4')
df = out_lipid_terms$Oli
summary1 = summary[rownames(df),]
summary1$AD = ifelse(summary1$niareagansc%in%c(1,2),1,0)

# unstratified
mod = model.matrix(~APOE4 + AD + age_death + msex + pmi, data=summary1)
fit <- lmFit(t(df[as.character(rownames(summary1)),]), design=mod)
fit <- eBayes(fit)
allgenesets <- topTable(fit, coef='AD', number=Inf, confint = T) %>% .[order(.$P.Value, decreasing = F),]
all_data[['res']][['lipid_associated']][['unstratified_AD_effect']] = allgenesets

# stratify by APOE status
summary2 = summary1[summary1$APOE=='e4',]
mod = model.matrix(~AD + age_death + msex + pmi, data=summary2)
fit <- lmFit(t(df[as.character(rownames(summary2)),]), design=mod)
fit <- eBayes(fit)
allgenesets <- topTable(fit, coef='AD', number=Inf, confint = T) %>% .[order(.$P.Value, decreasing = F),]
all_data[['res']][['lipid_associated']][['APOE4_stratified_AD_effect']] = allgenesets

# stratify by APOE status
summary3 = summary1[summary1$APOE=='e3',]
mod = model.matrix(~AD + age_death + msex + pmi, data=summary3)
fit <- lmFit(t(df[as.character(rownames(summary3)),]), design=mod)
fit <- eBayes(fit)
allgenesets <- topTable(fit, coef='AD', number=Inf, confint = T) %>% .[order(.$P.Value, decreasing = F),]
all_data[['res']][['lipid_associated']][['APOE33_stratified_AD_effect']] = allgenesets

# stratify by APOE status
summary4 = summary1[summary1$apoe_genotype==34,]
mod = model.matrix(~AD + age_death + msex + pmi, data=summary4)
fit <- lmFit(t(df[as.character(rownames(summary4)),]), design=mod)
fit <- eBayes(fit)
allgenesets <- topTable(fit, coef='AD', number=Inf, confint = T) %>% .[order(.$P.Value, decreasing = F),]
all_data[['res']][['lipid_associated']][['APOE34_stratified_AD_effect']] = allgenesets

# get all the effect sizes and p-values as a table
cholest_path = 'cholesterol biosynthesis III (via desmosterol) Homo sapiens PWY66-4'
out = list()
names = names(all_data[['res']][['lipid_associated']])
names = names[!names%in%c('gsva_out')]

for(i in names){
    print(i)
    f = all_data[['res']][['lipid_associated']][[i]][cholest_path,c('logFC', 'P.Value')]
    f$name = i
    out[[i]] = f
}

# save the table with logFC values and p-value
write.csv(do.call('rbind', out),'../data/other_analyses_outputs/stratified_stats.csv')
saveRDS(all_data, '../data/other_analyses_outputs/stratified_anaylsis.rds')
print('done')
