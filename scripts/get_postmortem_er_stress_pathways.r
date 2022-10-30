########## ER stress related to extended data figure 8 #############
##########################################################################
print('|| analysis of ER stress... ||')

source('../functions/pathway_setup.r')
source('../functions/pathway_fits.r')

##### required libraries ####
library('GSVA')
library("readxl")
library('SingleCellExperiment')
library('limma')
library('parallel')
library('tidyr')
library('stringr')

# load the pathways
pathways = readRDS('../data/other_analyses_outputs/pathways.rds')
all_paths = pathways$pathways$all

# output data
all_data = list()

# load the average expression data and metadata
expressed = readRDS('../data/single_cell_data/expressed_genes_per_celltype.rds')
av_expression = readRDS('../data/single_cell_data/individual_level_averages_per_celltype.rds')
summary = read.csv('../data/single_cell_data/metadata_by_individual.csv')
summary$APOE4 = ifelse(summary$apoe_genotype == '33', 0, 1)
rownames(summary) = summary[,'projid']

f1_data = readRDS('../data/other_analyses_outputs/pathway_scores.rds')
nebula = readRDS('../data/differentially_expressed_genes_data/E4_nebula_associations_by_celltype.rds')

# get ER stress/unfolded protein response pathways
paths = get_gset_names_by_category(c('unfolded protein'), names(all_paths))
paths = all_paths[paths]
low_removed = filter_lowly_exp_genes(expressed, paths)
order = c('Oli')
out = lapply(order, function(x) t(gsva(as.matrix(av_expression[[x]]), low_removed[[x]] , mx.diff=TRUE, verbose=TRUE, kcdf=c("Gaussian"), min.sz=5,parallel.sz=0)))
names(out) = order

# fit linear model to the gsva scores
predict = summary[as.character(rownames(out$Oli)),]
mod = model.matrix(~APOE4 + amyloid + nft + age_death + msex + pmi, data=predict)
fits = fit.gsva(mod, names(out), out, 'APOE4')

# for the pathway of interest, look at the specific genes
name = 'ATF6-mediated unfolded protein response (GO:0036500)'
df2 = as.data.frame(out$Oli[,name])
df2$genotype = summary[rownames(df2),'apoe_genotype']
df2$apoe = ifelse(df2$genotype == 33, 'E3', 'E4')
colnames(df2) = c('value', 'genotype', 'APOE')

df = nebula$Oli
df$padj = p.adjust(df$p_Apoe_e4yes, 'fdr')
degs = (df[df$padj<0.05,])
d = na.omit(degs[unname(unlist(low_removed$Oli[name])),])
x = d[,1]
names(x) = rownames(d)

all_data$ER_stress$ATF6_activity = df2
all_data$ER_stress$ATF6_pathway_degs = x
all_data$ER_stress$paths = low_removed

saveRDS(all_data, '../data/other_analyses_outputs/er_stress_results.rds')
print('done.')
