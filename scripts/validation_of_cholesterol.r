# load libraries
print('loading libraries..')
source('../functions/pathway_analyses.r')
library("readxl")
library('SingleCellExperiment')
library('GSVA')
library('limma')
library('parallel')
library('tidyr')
library('stringr')
library('ggplot2')
set.seed(5)

# load the nature data
nature_data = readRDS('../data/single_cell_data/Ctx.Mathys2019.aggregate.data.rds')
validation = list()

# load the pathway data
data = readRDS('../data/other_analyses_outputs/pathways.rds')
all_paths = data$pathways$all
f1_data = readRDS('../data/other_analyses_outputs/pathway_scores.rds')

# load the data
print('loading the data..')
expressed = readRDS('../data/single_cell_data/expressed_genes_per_celltype.rds')
av_expression = readRDS('../data/single_cell_data/individual_level_averages_per_celltype.rds')
summary = read.csv('../data/single_cell_data/metadata_by_individual.csv')
summary$APOE4 = ifelse(summary$apoe_genotype == '33','e3','e4')
rownames(summary) = summary[,'projid...2']

# only keep nature samples that are not overlapping with our dataset
meta = nature_data$metadata
meta_nonoverlapping = meta[!rownames(meta)%in%rownames(summary),]

# prepare the lipid pathways
validation_avs = nature_data$individual.level.averages

pathways = readRDS('../data/other_analyses_outputs/pathways.rds')
low_removed_lipid_paths = pathways$pathways$low_removed_lipid_associated

validation_cohort = meta_nonoverlapping[(!meta_nonoverlapping$apoe_genotype%in%c(23,24,0)),]

# filter by expressed genes (validation_avs)
exp = nature_data$detection.rate.by.celltype
expressed = lapply(colnames(exp), function(x) names(exp[,x][exp[,x]>0.05])) # using 5% because v2 --> lower detection rate
names(expressed) = colnames(exp)

# run GSVA
print('running GSVA...')
print('showing APOE4-associated pathways')
order = c('Oli')
out_go_terms = lapply(order, function(x) t(gsva(as.matrix(validation_avs[[x]][expressed[[x]],rownames(validation_cohort)]),all_paths[pathways$pathways$all_lipid_associated], mx.diff=TRUE, verbose=F, kcdf=c("Gaussian"), min.sz=5, max.sz = 500,parallel.sz=4)))
names(out_go_terms) = order
validation_cohort$APOE4 = ifelse(validation_cohort$apoe_genotype%in%c('24','34','44'),'e4','e3')
fits = get_fits(out_go_terms, validation_cohort)

validation[['APOE4_nature_fits']] = fits
df = as.data.frame(out_go_terms$Oli[,'Cholesterol metabolism'])
df$APOE = validation_cohort[rownames(df), 'APOE4']

colnames(df) = c('activity', 'APOE')
p <- ggplot(df, aes(x=as.factor(APOE), y=activity)) +
  geom_boxplot()
options(repr.plot.width = 2, repr.plot.height = 3, repr.plot.res = 300)

pdf('../plots/cholest_APOE4_nature_cohort.pdf', width = 2, height = 2)
p + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_classic()
dev.off()

# show the cholesterol genes
path = intersect(expressed[['Oli']], all_paths[['Cholesterol metabolism']])
write.csv(as.data.frame(path), '../data/supplementary_tables/cholesterol_metabolism_validation_apoe4_genes.csv')

# show validation data
x = c('nft', 'amyloid', 'pmi', 'age_death')

# show the metadata variables
print('plotting metadata variables...')
cols = c('grey', 'red')
names(cols) = c('33', '34')

for(i in x){
  pdf(paste0('../plots/nature_cohort_apoe_comparison_',i,'.pdf'), width = 2, height = 2)
  print(boxplot_w_stats(df = as.data.frame(validation_cohort), x = 'APOE4', y = i, palette = cols, comparisons = list(c('e3', 'e4')), xlab = '', ylab = i, width = .5, alpha = .5))
  dev.off()
}

# get fits for AD-associated pathways
print('showing AD-associated pathways')
df = read.csv('../data/single_cell_data/metadata_PFC_all_individuals_092520.tsv', sep = '\t')
rownames(df) = df$projid
validation_cohort$AD = ifelse(df[rownames(validation_cohort), 'niareagansc']%in%c(1,2) & validation_cohort$cogdx==4,1,ifelse(df[rownames(validation_cohort), 'niareagansc']%in%c(3,4) & validation_cohort$cogdx==1,0,0))
validation_cohort = validation_cohort[validation_cohort$APOE4=='e3',]
#validation_cohort = validation_cohort[validation_cohort$AD!='other',]
validation_cohort$AD = as.numeric(validation_cohort$AD)

mod = model.matrix(~AD + age_death  + msex + pmi, data=validation_cohort)
fit <- lmFit(t(out_go_terms$Oli[rownames(validation_cohort),]), design=mod)
fit <- eBayes(fit)
allgenesets_noAD_vs_AD <- topTable(fit, coef='AD', number=Inf, confint = T) %>% .[order(.$P.Value, decreasing = F),]

path = intersect(expressed[['Oli']], all_paths[['sterol biosynthetic process (GO:0016126)']])
validation[['AD_nature_fits']] = fits

df = as.data.frame(out_go_terms$Oli[rownames(validation_cohort),'sterol biosynthetic process (GO:0016126)'])
write.csv(as.data.frame(path), '../data/supplementary_tables/sterol_biosynthesis_validation_AD_genes.csv')
df$AD = validation_cohort[rownames(df), 'AD']

colnames(df) = c('activity', 'AD')
p <- ggplot(df, aes(x=as.factor(AD), y=activity)) +
  geom_boxplot()
options(repr.plot.width = 2, repr.plot.height = 3, repr.plot.res = 300)

pdf('../plots/cholest_APOE4_nature_cohort.pdf', width = 2, height = 2)
p + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_classic()
dev.off()

# show validation data
x = c('nft', 'amyloid', 'pmi', 'age_death')

# show the metadata variables
print('plotting metadata variables...')
cols = c('grey', 'red')
names(cols) = c('0', '1')

for(i in x){
  pdf(paste0('../plots/nature_cohort_apoe_comparison_',i,'.pdf'), width = 2, height = 2)
  print(boxplot_w_stats(df = as.data.frame(validation_cohort), x = 'AD', y = i, palette = cols, comparisons = list(c('0', '1')), xlab = '', ylab = i, width = .5, alpha = .5))
  dev.off()
}


saveRDS(validation, '../data/other_analyses_outputs/Nature_validation_fits.rds')
