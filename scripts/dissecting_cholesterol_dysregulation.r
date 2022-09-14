########## cholesterol analysis for figure 2 ##########
#######################################################
print('|| dissecting_cholesterol_dysregulation... ||')

source('../functions/pathway_setup.r')

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
low_removed_lipid_paths = pathways$pathways$low_removed_lipid_associated

# output data
all_data = list()

# load the average expression data and metadata
expressed = readRDS('../data/single_cell_data/expressed_genes_per_celltype.rds')
av_expression = readRDS('../data/single_cell_data/individual_level_averages_per_celltype.rds')
summary = read.csv('../data/single_cell_data/metadata_by_individual.csv')
summary$APOE4 = ifelse(summary$apoe_genotype == '33', 0, 1)
rownames(summary) = summary[,'projid...2']

f1_data = readRDS('../data/other_analyses_outputs/pathway_scores.rds')
nebula = readRDS('../data/differentially_expressed_genes_data/E4_nebula_associations_by_celltype.rds')$Oli

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
nebula$padj = p.adjust(nebula$p_Apoe_e4yes, 'fdr')

# get genes of interest
genes = unique(unname(unlist(paths)))
genes = genes[genes%in%expressed$Oli]
nebula_sele = nebula[genes,]
nebula_sele_sig = nebula_sele[nebula_sele$padj<0.05,]
x = nebula_sele_sig$logFC_Apoe_e4yes
names(x) = rownames(nebula_sele_sig)

all_data$deg_level_analysis$degs = x
all_data$deg_level_analysis$sterol_paths = paths

saveRDS(all_data, '../data/other_analyses_outputs/cholesterol_analysis.rds')
print('done')
