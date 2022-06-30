########## cholesterol analysis for figure 2 ##########
#######################################################

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

# get union of cholesterol biosynthesis genes and make union for density plot
print('getting union of cholesterol biosynthesis genes...')
oli.lipid.fits = f1_data$lipid$fits$Oli
lipid.paths.oli = low_removed_lipid_paths$Oli
paths = rownames(oli.lipid.fits[oli.lipid.fits$P.Value<0.05,])
biosynth_genes = unique(unname(unlist(lipid.paths.oli[paths])))
temp = list()
temp[['biosynth']] = biosynth_genes
out = gsva(as.matrix(av_expression[['Oli']]), temp, mx.diff=TRUE, verbose=F, kcdf=c("Gaussian"), min.sz=5, max.sz = 150, parallel.sz=50)
df = as.data.frame(t(out))
df$apoe_geno = summary[rownames(df), 'apoe_genotype']
df$AD = summary[rownames(df), 'niareagansc']
df$AD = ifelse(df$AD%in%c(1,2), 'AD', 'noAD')
df = as.data.frame(t(out))
df$apoe_geno = summary[rownames(df), 'apoe_genotype']
df$AD = summary[rownames(df), 'niareagansc']
df$AD = ifelse(df$AD%in%c(1,2), 'AD', 'noAD')

all_data[['union_cholest_biosynth']]$genes = biosynth_genes
all_data[['union_cholest_biosynth']]$density_plot_input = df

#cholesterol gene-level analysis
print('cholesterol gene-level analysis')
sterol_paths0 = get_gset_names_by_category(c('cholesterol storage', 'lipid storage', 'lipid droplets', 'sterol storage'), names(all_paths))
sterol_paths1 = get_gset_names_by_category(c('cholesterol transport','cholesterol export', 'lipid transport','lipid export','sterol transport', 'sterol export', 'cholesterol efflux', 'sterol efflux', 'lipid efflux'), names(all_paths))
sterol_paths2 = get_gset_names_by_category(c('cholesterol biosynthesis','cholesterol homeostasis', 'lipid biosynthesis','lipid homeostasis','sterol biosynthesis', 'sterol homeostasis'), names(all_paths))
paths = all_paths[unique(c(sterol_paths0,sterol_paths1,sterol_paths2))]

# get the degs
print('getting the degs..')
df = nebula$Oli
df$padj = p.adjust(df$p_Apoe_e4yes, 'fdr')
degs = rownames(df[df$padj<0.05,])
degs = degs[degs%in%expressed$Oli]

# get the hmap input
print('getting the hmap input...')
outs = list()
sets = paths
for(i in names(sets)){
    curr = as.data.frame(sets[[i]])
    curr$path = i
    curr$present = 1
    outs[[i]] = curr
}
outs = do.call('rbind',outs)
outs = outs[outs[,1]!='',]

d = as.data.frame(pivot_wider(outs, values_from = 'present', names_from = 'path'))
rownames(d) = d[,1]
d[,1] = NULL
plot = d[rownames(d)%in%degs,]
plot[is.na(plot)] = 0
plot = plot[,colSums(plot)>0]

# show as barplot
dff = df[rownames(plot),1,drop = F]
x = dff[,1]
names(x) = rownames(dff)

all_data$deg_level_analysis$degs = x
all_data$deg_level_analysis$sterol_paths = paths

saveRDS(all_data, '../data/other_analyses_outputs/cholesterol_analysis.rds')
print('done')
