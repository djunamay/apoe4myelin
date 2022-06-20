##### required libraries ####
source('../functions/pathway_analyses.r')
library('GSVA')
library("readxl")
library('SingleCellExperiment')
library('limma')
library('parallel')
library('tidyr')
library('stringr')
library('pbmcapply')
set.seed(5)

# load observed pathway scores
scores = readRDS('../data_outputs/pathway_scores.rds')
df = (scores$res$apoe_associated$scores$all$Oli)
observation = nrow(df[df$P.Value<0.05,])

# load the pathways
print('loading pathways')
pathways = readRDS('../data_outputs/pathways.rds')
low_removed_bp = pathways$pathways$low_removed_bp
apoe_gsets_low_removed = pathways$pathways$apoe_gsets_low_removed
apoe_gsets = pathways$pathways$apoe_gsets_all
all_paths = pathways$pathways$all
low_removed_lipid_paths = pathways$pathways$low_removed_lipid_associated

# load the average expression data and metadata
print('loading data') # TODO simplify data
e = readRDS('../../submission_code_09012021/data/Summary.DE.celltype.rds')
expressed = metadata(e)$expressed.genes
av_expression = readRDS('../../submission_code_09012021/data/Averages.by.celltype.by.individual.rds')

summary = as.data.frame(read_excel('../../submission_code_09012021/data/Overview_PFC_projids_data_384.xlsx'))
summary = summary[!duplicated(summary[,'projid...2']),]
rownames(summary) = summary[['projid...2']]

# get scaffolding for gset sampling

summary$APOE4 = ifelse(summary$apoe_genotype == '33','e3','e4')


get_rand_res_permute_gsets = function(lengths_all, all_genes_all, celltype, av_expression, summary){
    lengths = lengths_all[[celltype]]
    temp = all_genes_all[[i]]
    rand_gsets = lapply(1:length(lengths), function(x) sample(temp, replace = F, lengths[x]))
    out_lipid_terms = lapply(celltype, function(x) t(gsva(as.matrix(av_expression[[x]]), rand_gsets, mx.diff=TRUE, verbose=F, kcdf=c("Gaussian"), min.sz=5, max.sz = 150, parallel.sz = 0)))
    names(out_lipid_terms) = celltype
    fits = get_fits(out_lipid_terms, summary)
    temp = (fits[[celltype]])
    rand = nrow(temp[temp$P.Value<0.05,])
    names(rand) = celltype
    return(rand)
}

# get the null results
print('simulating null model')
order = c('Oli', 'Opc', 'Mic', 'Ast', 'Ex', 'In')
print('starts')


lengths_all = list()
all_genes_all = list()
for(i in order){
    paths = low_removed_lipid_paths[[i]]
    all_genes_all[[i]] = unique(unname(unlist(paths)))
    lengths_all[[i]] = unlist(lapply(1:length(paths), function(x) length(paths[[x]])))

}

x = unlist(pbmclapply(1:10000, function(y) lapply(order, function(x) get_rand_res_permute_gsets(lengths_all, all_genes_all, x, av_expression, summary)), mc.cores = 40))
saveRDS(x, '../data_outputs/gset_permutation_results.rds')
