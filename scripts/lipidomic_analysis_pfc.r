library(readxl)

# load the che data
data = as.data.frame(read_excel('../data/lipidomics_dataset/pfc_lipidomics/ChE summary_cyc_05312022_all samples.xlsx'))
data_new = na.omit(data[(as.numeric(data[,'Ratio of signal to noise'])>3 & as.numeric(data[,'Peak Quality'])>=0.6),])
index = endsWith(data_new$label,')') | endsWith(data_new$label,'1')
data_new = data_new[index,]
data_subset = data_new[,9:ncol(data_new)]
rownames(data_subset) = data_new$label
data_new = data_new[data_new$label!='ChE(18:1D7)_1',]
data_subset = data_subset[rownames(data_subset)!='ChE(18:1D7)_1',]

# load the metadata
biospecimen = read.csv('../data/lipidomics_dataset/pfc_lipidomics/ROSMAP_Lipidomics_Emory_biospecimen_metadata.csv')
rownames(biospecimen) = biospecimen$specimenID
full_meta = read.csv('../data/lipidomics_dataset/pfc_lipidomics/ROSMAP_clinical.csv')
rownames(full_meta) = full_meta$individualID

data_subset = data_subset[,biospecimen[colnames(data_subset),'tissue']=='dorsolateral prefrontal cortex']
colnames(data_subset) = full_meta[biospecimen[colnames(data_subset),'individualID'],'projid']
data_subset = as.data.frame(t(data_subset))
rownames(full_meta) = full_meta$projid
data_subset$apoe_genotype = full_meta[rownames(data_subset), 'apoe_genotype']

nia = read.csv('../data/lipidomics_dataset/pfc_lipidomics/metadata_PFC_all_individuals_092520.tsv', sep = '\t')
full_meta = merge(full_meta, nia, by = 'projid')
rownames(full_meta) = full_meta$projid
data_subset$niareagansc = full_meta[rownames(data_subset), 'niareagansc']
data_subset$AD = ifelse(data_subset$niareagansc%in%c(1,2), 'AD', 'noAD')
data_subset = data_subset[data_subset$apoe_genotype%in%c(33,34,44),]

# perform wilcoxon test for each lipid species of interest
stats = list()
for(i in colnames(data_subset)[1:3]){
  y = data_subset[data_subset$apoe_genotype=='33',i]
  x = data_subset[data_subset$apoe_genotype!='33',i]
  temp = wilcox.test(x,y,conf.int=T, conf.level = 0.90)
  stats[[i]] = c(temp$conf.int, temp$estimate, temp$p.value)
}

stats = do.call('rbind', stats)
stats = as.data.frame(stats)
colnames(stats) = c('lower_90_CI', 'upper_90_CI', 'sample_estimates(median_difference_44_vs_33)', 'p.value')
stats$grp = 'unstratified'

stats_noAD = list()
data_subset_noAD = data_subset[data_subset$AD=='noAD',]
for(i in colnames(data_subset_noAD)[1:3]){
  y = data_subset_noAD[data_subset_noAD$apoe_genotype=='33',i]
  x = data_subset_noAD[data_subset_noAD$apoe_genotype!='33',i]
  temp = wilcox.test(x,y,conf.int=T, conf.level = 0.90)
  stats_noAD[[i]] = c(temp$conf.int, temp$estimate, temp$p.value)
}

stats_noAD = do.call('rbind', stats_noAD)
stats_noAD = as.data.frame(stats_noAD)
colnames(stats_noAD) = c('lower_90_CI', 'upper_90_CI', 'sample_estimates(median_difference_44_vs_33)', 'p.value')
stats_noAD$grp = 'noAD'

stats_AD = list()
data_subset_AD = data_subset[data_subset$AD=='AD',]
for(i in colnames(data_subset_AD)[1:3]){
  y = data_subset_AD[data_subset_AD$apoe_genotype=='33',i]
  x = data_subset_AD[data_subset_AD$apoe_genotype!='33',i]
  temp = wilcox.test(x,y,conf.int=T, conf.level = 0.90)
  stats_AD[[i]] = c(temp$conf.int, temp$estimate, temp$p.value)
}

stats_AD = do.call('rbind', stats_AD)
stats_AD = as.data.frame(stats_AD)
colnames(stats_AD) = c('lower_90_CI', 'upper_90_CI', 'sample_estimates(median_difference_44_vs_33)', 'p.value')
stats_AD$grp = 'AD'

df = rbind(stats, stats_noAD, stats_AD)
write.csv(df, '../data/other_analyses_outputs/cc_lipidomics_data.csv')

# also save the metadata --> add age etc to it
write.csv(df, '../data/other_analyses_outputs/cc_lipidomics_data.csv')

# make metadata distribution plots

# show the plots
