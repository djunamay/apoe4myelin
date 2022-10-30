########## lipidomic analysis for extended data figure 4 #############
#######################################################################
print('|| pfc_lipidomic_analysisis - get data... ||')

source('../functions/differential_expression.r')

# required packages
library(readxl)
library(ggplot2)
library(ggpubr)

# load the che data
data = as.data.frame(read_excel('../data/pfc_lipidomics/ChE summary_cyc_05312022_all samples.xlsx'))
data_new = na.omit(data[(as.numeric(data[,'Ratio of signal to noise'])>3 & as.numeric(data[,'Peak Quality'])>=0.6),]) # filter species
index = endsWith(data_new$label,')') | endsWith(data_new$label,'1') | endsWith(data_new$label,'2') # filter species (don't consider peak combinations _c)
data_new = data_new[index,] # filter species
data_subset = data_new[,9:ncol(data_new)] # get only the lipidomic values
rownames(data_subset) = data_new$label
data_subset = data_subset[!rownames(data_subset)%in%c('ChE(18:1D7)_2','ChE(18:1D7)_1'),] # remove spike ins

# load the metadata
biospecimen = read.csv('../data/pfc_lipidomics/ROSMAP_Lipidomics_Emory_biospecimen_metadata.csv')
rownames(biospecimen) = biospecimen$specimenID
full_meta = read.csv('../data/pfc_lipidomics/ROSMAP_clinical.csv')
rownames(full_meta) = full_meta$individualID
name_conversion = merge(full_meta, biospecimen, by = 'individualID')[,c('projid', 'specimenID', 'individualID')]
rownames(name_conversion) = name_conversion$specimenID

# renaming
data_subset = data_subset[,biospecimen[colnames(data_subset),'tissue']=='dorsolateral prefrontal cortex'] # only keep PFC data
data_subset = as.data.frame(t(data_subset))
data_subset = merge(data_subset, name_conversion, by = 0)

# add additional variable info to the data
nia = read.csv('../data/pfc_lipidomics/metadata_PFC_all_individuals_092520.tsv', sep = '\t')
full_meta = merge(full_meta, nia, by = 'projid')
rownames(full_meta) = full_meta$projid

data_subset$niareagansc = full_meta[as.character(data_subset$projid), 'niareagansc']
data_subset$pmi = full_meta[as.character(data_subset$projid), 'pmi.x']
data_subset$nft = full_meta[as.character(data_subset$projid), 'nft']
data_subset$amyloid = full_meta[as.character(data_subset$projid), 'amyloid']
data_subset$age_death = full_meta[as.character(data_subset$projid), 'age_death.y']
data_subset$msex = full_meta[as.character(data_subset$projid), 'msex.x']
data_subset$apoe_genotype = full_meta[as.character(data_subset$projid), 'apoe_genotype.x']
data_subset$AD = ifelse(data_subset$niareagansc%in%c(1,2), 'AD', 'noAD')
data_subset = na.omit(data_subset)
data_subset = data_subset[data_subset$apoe_genotype%in%c(33,34,44),] # only consider 33 and 34 and 44 individuals

# perform wilcoxon test for each lipid species of interest
stats = lipid_species_wilcox_test(data_subset,'unstratified')
stats_noAD = lipid_species_wilcox_test(data_subset[data_subset$niareagansc%in%c('3','4'),],'noAD')
stats_AD = lipid_species_wilcox_test(data_subset[data_subset$niareagansc%in%c('1','2'),],'AD')

df = rbind(stats, stats_noAD, stats_AD)
df$species = rep(c( "ChE(18:1)_1", "ChE(18:1)_2", "ChE(18:2)_1", "ChE(20:4)_1"),3)
write.csv(df, '../data/supplementary_tables/pfc_lipidomics_data.csv')

# also save the metadata
write.csv(data_subset, '../data/supplementary_tables/pfc_lipidomics_data_metadata.csv')

# add the extra supp tables
stats = data[,1:8]
write.csv(stats, '../data/supplementary_tables/pfc_lipidomics_qc_metrics.csv')

all_data_no_qc = data[,9:ncol(data)]
rownames(all_data_no_qc) = stats$label
write.csv(all_data_no_qc, '../data/supplementary_tables/lipicomics_all_data_no_qc.csv')

print('done.')
