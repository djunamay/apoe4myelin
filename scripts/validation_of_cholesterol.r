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
nature_data = readRDS('../data/Ctx.Mathys2019.aggregate.data.rds')

validation = list()
# load the pathway data
data = readRDS('../data/pathways.rds')
all_paths = data$pathways$all
f1_data = readRDS('../data_outputs/pathway_scores.rds')
#f2_data = readRDS('../data/figure_2_data.rds')

# load all the metadata
summary_train = as.data.frame(read_excel('../../submission_code_09012021/data/Overview_PFC_projids_data_384.xlsx'))
summary_train = summary_train[!duplicated(summary_train[,'projid...2']),]
rownames(summary_train) = summary_train[['projid...2']]
av_expression = readRDS('../../submission_code_09012021/data/Averages.by.celltype.by.individual.rds')
summary_train2 = summary_train[colnames(av_expression$Oli),]

# only keep nature samples that are not overlapping with our dataset
meta = nature_data$metadata
meta_nonoverlapping = meta[!rownames(meta)%in%rownames(summary_train2),]

# prepare the lipid pathways
validation_avs = nature_data$individual.level.averages
exp = nature_data$detection.rate.by.celltype

# test if the same lipid pathways are dysregulated
expressed = lapply(colnames(exp), function(x) names(exp[,x][exp[,x]>0.05])) # using 5% because v2 --> lower detection rate
names(expressed) = colnames(exp)

pathways = readRDS('../data_outputs/pathways.rds')
low_removed_lipid_paths = pathways$pathways$low_removed_lipid_associated

########## Testing E4 effect on lipid pathways
# get validation metadata
validation_cohort = meta_nonoverlapping[(!meta_nonoverlapping$apoe_genotype%in%c(23,0)),]
validation_cohort = summary_train[as.character(rownames(validation_cohort)),]
validation_cohort = validation_cohort[validation_cohort$niareagansc%in%c(1,2),]

# look at the metadata variable distributions - get the key metadata dist
box::use(../../../ABCA7_LOF_2022/scmod_R/general_functions)

# pick the cohort stratification that has the most evenly matched covariates
cols = c('grey', 'red')
names(cols) = c('33', '34')

# show validation data
x = c('nft', 'amyloid', 'pmi', 'age_death')

for(i in x){
  pdf(paste0('../plots/nature_cohort_apoe_comparison_',i,'.pdf'), width = 2, height = 2)
  print(general_functions$boxplot_w_stats(df = as.data.frame(validation_cohort), x = 'apoe_genotype', y = i, palette = cols, comparisons = list(c('33', '34')), xlab = '', ylab = i, width = .5, alpha = .5))
  dev.off()
}

# print the number of males per genotype
print('msex by apoe genotype = ')
print(table(validation_cohort$msex, validation_cohort$apoe_genotype))

print('niareagansc by apoe genotype = ')
print(table(validation_cohort$niareagansc, validation_cohort$apoe_genotype))

# run GSVA on reactome pathways
print('running GSVA...')
order = c('Oli')
out_go_terms = lapply(order, function(x) t(gsva(as.matrix(validation_avs[[x]][,rownames(validation_cohort)]),low_removed_lipid_paths[[x]], mx.diff=TRUE, verbose=F, kcdf=c("Gaussian"), min.sz=5, max.sz = 500,parallel.sz=4)))
names(out_go_terms) = order

# get linear model fits
print('getting linear model fits...')
out_go_terms_ordered = lapply(names(out_go_terms), function(x) out_go_terms[[x]][rownames(validation_cohort),])
names(out_go_terms_ordered) = names(out_go_terms)

validation_cohort$APOE4 = ifelse(validation_cohort$apoe_genotype%in%c('23','33'),'e3','e4')
fits = get_fits(out_go_terms_ordered, validation_cohort)
validation[['APOE4_nature_fits']] = fits
df = as.data.frame(out_go_terms$Oli[,'Regulation of cholesterol biosynthesis by SREBP (SREBF) Homo sapiens R-HSA-1655829'])
df$APOE = validation_cohort[rownames(df), 'APOE4']

colnames(df) = c('activity', 'APOE')
p <- ggplot(df, aes(x=as.factor(APOE), y=activity)) +
  geom_boxplot()
options(repr.plot.width = 2, repr.plot.height = 3, repr.plot.res = 300)

pdf('../plots/cholest_APOE4_nature_cohort.pdf', width = 2, height = 2)
p + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_classic()
dev.off()

########## Testing AD effect on lipid pathways

# do the same for AD
validation_cohort = meta_nonoverlapping[(meta_nonoverlapping$apoe_genotype%in%c(33)),]
validation_cohort = summary_train[as.character(rownames(validation_cohort)),]
validation_cohort = validation_cohort[(validation_cohort$niareagansc%in%c(1,2) & validation_cohort$cogdx%in%c(4,5)) | (validation_cohort$niareagansc%in%c(3,4) & validation_cohort$cogdx%in%c(1,2)),]
validation_cohort$AD = ifelse(validation_cohort$niareagansc%in%c(1,2) & validation_cohort$cogdx%in%c(4,5),1,0)

# run GSVA on reactome pathways
print('running GSVA...')
order = c('Oli')
lipid_paths_low_removed = f1_data$pathways$lipid_associated
out_lipid_terms = lapply(order, function(x) t(gsva(as.matrix(validation_avs[[x]]), low_removed_lipid_paths[[x]], mx.diff=TRUE, verbose=F, kcdf=c("Gaussian"), min.sz=5, max.sz = 1000)))
names(out_lipid_terms) = order
#low_removed_lipid_paths[[x]]

# get linear model fits
print('getting linear model fits...')
out_lipid_terms$Oli = out_lipid_terms$Oli[rownames(validation_cohort),]

validation_cohort$AD = ifelse(validation_cohort$niareagansc%in%c(1,2) & validation_cohort$cogdx%in%c(4,5),1,0)
mod = model.matrix(~AD + age_death  + msex + pmi, data=validation_cohort)
fit <- lmFit(t(out_lipid_terms$Oli), design=mod)
fit <- eBayes(fit)
allgenesets_noAD_vs_AD <- topTable(fit, coef='AD', number=Inf, confint = T) %>% .[order(.$P.Value, decreasing = F),]
validation[['AD_nature_fits']] = as.data.frame(allgenesets_noAD_vs_AD)
df = as.data.frame(out_lipid_terms$Oli[,'cholesterol biosynthesis III (via desmosterol) Homo sapiens PWY66-4'])
#df$APOE = validation_cohort[rownames(df), 'APOE4']
df$AD = validation_cohort[rownames(df), 'AD']

colnames(df) = c('activity', 'AD')
p <- ggplot(df, aes(x=as.factor(AD), y=activity)) +
  geom_boxplot()
options(repr.plot.width = 2, repr.plot.height = 3, repr.plot.res = 300)

pdf('../plots/cholest_AD_nature_cohort.pdf', width = 2, height = 2)
p + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_classic()
dev.off()

x = c('nft', 'amyloid', 'pmi', 'age_death')
cols = c('grey', 'red')
names(cols) = c('0', '1')

for(i in x){
    pdf(paste0('../plots/nature_cohort_AD_comparison_',i, '.pdf'), width = 2, height = 2)
    print(general_functions$boxplot_w_stats(df = as.data.frame(validation_cohort), x = 'AD', y = i, palette = cols, comparisons = list(c('0', '1')), xlab = '', ylab = i, width = .5, alpha = .5))
    dev.off()
}

# print the number of males per genotype
print('msex by apoe genotype = ')
print(table(validation_cohort$msex, validation_cohort$niareagansc))

print('niareagansc by apoe genotype = ')
print(table(validation_cohort$apoe_genotype, validation_cohort$niareagansc))

saveRDS(validation, '../data_outputs/validation_out.rds')
print('done')
