
source('../functions/pathway_analyses.r')

##### required libraries ####
library('GSVA')
library("readxl")
library('SingleCellExperiment')
library('limma')
library('parallel')
library('tidyr')
library('stringr')
set.seed(5)

# load the pathways
print('loading pathways')
pathways = readRDS('../data/other_analyses_outputs/pathways.rds')
low_removed_bp = pathways$pathways$low_removed_bp
apoe_gsets_low_removed = pathways$pathways$apoe_gsets_low_removed
apoe_gsets = pathways$pathways$apoe_gsets_all
all_paths = pathways$pathways$all
low_removed_lipid_paths = pathways$pathways$low_removed_lipid_associated

# output data
all_data = list()

# load the average expression data and metadata
print('loading data')
print('loading the data..')
expressed = readRDS('../data/single_cell_data/expressed_genes_per_celltype.rds')
av_expression = readRDS('../data/single_cell_data/individual_level_averages_per_celltype.rds')
summary = read.csv('../data/single_cell_data/metadata_by_individual.csv')
summary$APOE4 = ifelse(summary$apoe_genotype == '33','e3','e4')
rownames(summary) = summary[,'projid...2']

f1_data = readRDS('../data/other_analyses_outputs/pathway_scores.rds')
nebula = readRDS('../data/differentially_expressed_genes_data/E4_nebula_associations_by_celltype.rds')

#########################################
####### ER pathway-level analysis #######
#########################################
print('ER pathway-level analysis')

paths = get_gset_names_by_category(c('unfolded protein'), names(all_paths))
paths = all_paths[paths]
low_removed = filter_lowly_exp_genes(expressed, paths)
order = c('Oli')
out = lapply(order, function(x) t(gsva(as.matrix(av_expression[[x]]), low_removed[[x]] , mx.diff=TRUE, verbose=TRUE, kcdf=c("Gaussian"), min.sz=5,parallel.sz=0)))
names(out) = order
# check that gsva output rownames are identical
check_rownames(out)

# fit linear model to the gsva scores
print('fitting linear model...')
summary$APOE4 = ifelse(summary$apoe_genotype == '33',0,1)

predict = summary[as.character(rownames(out$Oli)),]
mod = model.matrix(~APOE4 + amyloid + nft + age_death + msex + pmi, data=predict)

fits = fit.gsva(mod, names(out), out, 'APOE4')

# look at the ER pathway more specifically
name = 'ATF6-mediated unfolded protein response (GO:0036500)'

# show the modules as boxplots
print('showing pathway as boxplot...')
df2 = as.data.frame(out$Oli[,name])
df2$genotype = summary[rownames(df2),'apoe_genotype']

df2$apoe = ifelse(df2$genotype == 33, 'E3', 'E4')
colnames(df2) = c('value', 'genotype', 'APOE')

print('getting degs...')
df = nebula$Oli
df$padj = p.adjust(df$p_Apoe_e4yes, 'fdr')
degs = (df[df$padj<0.05,])
browser()
d = na.omit(degs[unname(unlist(low_removed$Oli[name])),])
x = d[,1]
names(x) = rownames(d)

all_data$ER_stress$ATF6_activity = df2
all_data$ER_stress$ATF6_pathway_degs = x
all_data$ER_stress$paths = low_removed

saveRDS(all_data, '../data/other_analyses_outputs/er_stress_results.rds')
print('done.')
