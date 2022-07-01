########## lipidomic analysis for extended data figure 4 #############
##########################################################################
source('../functions/pathway_analyses.r')

# required packages
library(readxl)
library(ggplot2)
library(ggpubr)

# load the che data
data = as.data.frame(read_excel('../data/lipidomics_dataset/pfc_lipidomics/ChE summary_cyc_05312022_all samples.xlsx'))
data_new = na.omit(data[(as.numeric(data[,'Ratio of signal to noise'])>3 & as.numeric(data[,'Peak Quality'])>=0.5),]) # filter species
index = endsWith(data_new$label,')') | endsWith(data_new$label,'1') | endsWith(data_new$label,'2') # filter species
data_new = data_new[index,] # filter species
data_subset = data_new[,9:ncol(data_new)] # get only the lipidomic values
rownames(data_subset) = data_new$label
data_subset = data_subset[!rownames(data_subset)%in%c('ChE(18:1D7)_2','ChE(18:1D7)_1'),] # remove spike ins

# load the metadata
biospecimen = read.csv('../data/lipidomics_dataset/pfc_lipidomics/ROSMAP_Lipidomics_Emory_biospecimen_metadata.csv')
rownames(biospecimen) = biospecimen$specimenID
full_meta = read.csv('../data/lipidomics_dataset/pfc_lipidomics/ROSMAP_clinical.csv')
rownames(full_meta) = full_meta$individualID

# renaming
data_subset = data_subset[,biospecimen[colnames(data_subset),'tissue']=='dorsolateral prefrontal cortex'] # only keep PFC data
biospecimen_id = colnames(data_subset) # save the biospecimen IDs
colnames(data_subset) = full_meta[biospecimen[colnames(data_subset),'individualID'],'projid'] # rename biospecimens according to projid
data_subset = as.data.frame(t(data_subset))
data_subset$biospecimen = biospecimen_id # add back the biospecimen ID

# add apoe genotype to data and subset the metadata
rownames(full_meta) = full_meta$projid
data_subset$apoe_genotype = full_meta[rownames(data_subset), 'apoe_genotype']
full_meta = full_meta[rownames(data_subset),]

# add additional variable info to the data
nia = read.csv('../data/lipidomics_dataset/pfc_lipidomics/metadata_PFC_all_individuals_092520.tsv', sep = '\t')
full_meta = merge(full_meta, nia, by = 'projid')
rownames(full_meta) = full_meta$projid
data_subset$niareagansc = full_meta[rownames(data_subset), 'niareagansc']
data_subset$pmi = full_meta[rownames(data_subset), 'pmi.x']
data_subset$nft = full_meta[rownames(data_subset), 'nft']
data_subset$amyloid = full_meta[rownames(data_subset), 'amyloid']
data_subset$age_death = full_meta[rownames(data_subset), 'age_death.y']
data_subset$msex = full_meta[rownames(data_subset), 'msex.y']
data_subset = na.omit(data_subset)
data_subset$AD = ifelse(data_subset$niareagansc%in%c(1,2), 'AD', 'noAD')

# perform wilcoxon test for each lipid species of interest
stats = list()
for(i in colnames(data_subset)[1:3]){
  y = data_subset[data_subset$apoe_genotype%in%c('33','23','22'),i]
  x = data_subset[!data_subset$apoe_genotype%in%c('33','23','22'),i]
  temp = wilcox.test((x),(y),conf.int=T, conf.level = 0.90)
  stats[[i]] = c(temp$conf.int, temp$estimate, temp$p.value)
}

stats = do.call('rbind', stats)
stats = as.data.frame(stats)
colnames(stats) = c('lower_90_CI', 'upper_90_CI', 'sample_estimates.median_difference_APOE4carrier_vs_noncarrier.', 'p.value')
stats$grp = 'unstratified'

stats_noAD = list()
data_subset_noAD = data_subset[data_subset$niareagansc=='3',]
for(i in colnames(data_subset_noAD)[1:3]){
  y = data_subset_noAD[data_subset_noAD$apoe_genotype%in%c('33','23','22'),i]
  x = data_subset_noAD[!data_subset_noAD$apoe_genotype%in%c('33','23','22'),i]
  temp = wilcox.test((x),(y),conf.int=T, conf.level = 0.90)
  stats_noAD[[i]] = c(temp$conf.int, temp$estimate, temp$p.value)
}

stats_noAD = do.call('rbind', stats_noAD)
stats_noAD = as.data.frame(stats_noAD)
colnames(stats_noAD) = c('lower_90_CI', 'upper_90_CI', 'sample_estimates.median_difference_APOE4carrier_vs_noncarrier.', 'p.value')
stats_noAD$grp = 'noAD'

stats_AD = list()
data_subset_AD = data_subset[data_subset$niareagansc%in%c('1','2'),]
for(i in colnames(data_subset_AD)[1:3]){
  y = data_subset_AD[data_subset_AD$apoe_genotype%in%c('33','23','22'),i]
  x = data_subset_AD[!data_subset_AD$apoe_genotype%in%c('33','23','22'),i]
  temp = wilcox.test((x),(y),conf.int=T, conf.level = 0.90)
  stats_AD[[i]] = c(temp$conf.int, temp$estimate, temp$p.value)
}

stats_AD = do.call('rbind', stats_AD)
stats_AD = as.data.frame(stats_AD)
colnames(stats_AD) = c('lower_90_CI', 'upper_90_CI', 'sample_estimates.median_difference_APOE4carrier_vs_noncarrier.', 'p.value')
stats_AD$grp = 'AD'

df = rbind(stats, stats_noAD, stats_AD)
df$species = c( "ChE(16:1)" ,"ChE(18:1)_1" , "ChE(18:1)_2",  "ChE(16:1)" ,  "ChE(18:1)_1",
"ChE(18:1)_2" ,"ChE(16:1)" ,  "ChE(18:1)_1" ,"ChE(18:1)_2")
write.csv(df, '../data/other_analyses_outputs/pfc_lipidomics_data.csv')

# also save the metadata
write.csv(data_subset, '../data/other_analyses_outputs/pfc_lipidomics_data_metadata.csv')

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

# add the stacked bar plot
pdf('../plots/Extended_4/amyloid_msex_lipid_cohort.pdf')
get_barplot(data_subset, 'grp', 'msex')
dev.off()

# add the extra supp tables
stats = data[,1:8]
write.csv(stats, '../data/other_analyses_outputs/pfc_lipidomics_qc_metrics.csv')

all_data_no_qc = data[,9:ncol(data)]
rownames(all_data_no_qc) = stats$label
write.csv(all_data_no_qc, '../data/other_analyses_outputs/lipicomics_all_data_no_qc.csv')


# show the plots for lipid species of interest (pval < 0.1)
df[df$p.value<0.1,]
data_subset$APOE = factor(data_subset$APOE, levels = c('APOE4-noncarrier','APOE4-carrier'))
data_subset$AD = factor(data_subset$AD, levels = c('noAD', 'AD'))
pdf('../plots/Extended_4/ChE(16:1).pdf', width = 10, height = 5)
ggplot(data_subset, aes(x=(APOE), y=(data_subset[,'ChE(16:1)']), col = APOE)) +
stat_compare_means(comparisons = list(c('APOE4-noncarrier', 'APOE4-carrier')),method = "wilcox.test") + geom_violin() +  geom_boxplot(width = .1) + geom_jitter(size = .5) + theme_classic() + facet_grid(~ AD)
dev.off()

pdf('../plots/Extended_4/ChE(18:1)_2.pdf', width = 10, height = 5)
ggplot(data_subset, aes(x=(APOE), y=(data_subset[,'ChE(18:1)_2']), col = APOE)) +
stat_compare_means(comparisons = list(c('APOE4-noncarrier', 'APOE4-carrier')),method = "wilcox.test") + geom_violin() +  geom_boxplot(width = .1) + geom_jitter(size = .5) + theme_classic() + facet_grid(~ AD)
dev.off()

pdf('../plots/che_16_1_boxplot.pdf', width = 4, height = 7)
ggplot(data_subset, aes(x=(APOE), y=data_subset[,'ChE(16:1)'], col = APOE)) +
stat_compare_means(comparisons = list(c('APOE4-noncarrier', 'APOE4-carrier')),method = "wilcox.test") +  geom_boxplot(width = 1) + theme_classic() + facet_grid(~ AD)
dev.off()

pdf('../plots/che_18_2_boxplot.pdf', width = 4, height = 7)
ggplot(data_subset, aes(x=(APOE), y=data_subset[,'ChE(18:1)_2'], col = APOE)) +
stat_compare_means(comparisons = list(c('APOE4-noncarrier', 'APOE4-carrier')),method = "wilcox.test") +  geom_boxplot(width = 1) + theme_classic() + facet_grid(~ AD)
dev.off()

print('done.')
