########## script 4 in APOE4 impairs myelination via cholesterol dysregulation in oligodendrocytes ##########
#############################################################################################################

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
pathways = readRDS('../data/pathways.rds')
low_removed_bp = pathways$pathways$low_removed_bp
apoe_gsets_low_removed = pathways$pathways$apoe_gsets_low_removed
apoe_gsets = pathways$pathways$apoe_gsets_all
all_paths = pathways$pathways$all
low_removed_lipid_paths = pathways$pathways$low_removed_lipid_associated

# output data
all_data = list()

# load the average expression data and metadata
print('loading data') # TODO simplify data
e = readRDS('../../submission_code_09012021/data/Summary.DE.celltype.rds')
expressed = metadata(e)$expressed.genes
av_expression = readRDS('../../submission_code_09012021/data/Averages.by.celltype.by.individual.rds')

summary = as.data.frame(read_excel('../../submission_code_09012021/data/Overview_PFC_projids_data_384.xlsx'))
summary = summary[!duplicated(summary[,'projid...2']),]
rownames(summary) = summary[['projid...2']]

f1_data = readRDS('../data/pathway_scores.rds')

nebula = readRDS('../../submission_code_09012021/data/E4_nebula_associations_by_celltype.rds')

###########################################################################################
####### get union of cholesterol biosynthesis genes and make union for density plot #######
###########################################################################################
print('getting union of cholesterol biosynthesis genes...')
oli.lipid.fits = f1_data$res$lipid_associated$fits$Oli
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

###############################################
####### cholesterol gene-level analysis #######
###############################################
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

saveRDS(all_data, '../data/cholesterol_analysis.rds')
print('done')
