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
library(ggplot2)
set.seed(5)

# load observed pathway scores
scores = readRDS('../data/other_analyses_outputs/pathway_scores.rds')
paths = list()
medians = list()
for(i in c('Oli', 'Opc', 'Ex', 'In', 'Ast', 'Mic')){
  df = (scores$APOE$scores_all$all[[i]])
  names = df[df$P.Value<0.05,'names']
  medians[[i]] = median((df[df$P.Value<0.05,'logFC']))
  paths[[i]] = names
}

# load the pathways
print('loading pathways')
pathways = readRDS('../data/other_analyses_outputs/pathways.rds')
apoe_gsets_low_removed = pathways$pathways$apoe_gsets_low_removed

# load the average expression data and metadata
print('loading the data..')
expressed = readRDS('../data/single_cell_data/expressed_genes_per_celltype.rds')
av_expression = readRDS('../data/single_cell_data/individual_level_averages_per_celltype.rds')
summary = read.csv('../data/single_cell_data/metadata_by_individual.csv')
summary$APOE4 = ifelse(summary$apoe_genotype == '33','e3','e4')
rownames(summary) = summary[,'projid...2']

all_data = readRDS('../data/other_analyses_outputs/pathway_scores.rds')
out_apoe_terms = all_data$APOE$gsva_out
summary$APOE4 = ifelse(summary$apoe_genotype == '33','e3','e4')

get_rand_res_permute_labels = function(out_lipid_terms, summary){
    summary$APOE4 = sample(summary$APOE4)
    fits = get_fits(out_lipid_terms, summary)
    rand = lapply(names(fits), function(x) median(fits[[x]][paths[[x]],'logFC']))
    names(rand) = names(fits)
    return(unlist(rand))
}

out = pbmclapply(1:10000, function(x) get_rand_res_permute_labels(out_apoe_terms, summary), mc.cores = 16)

# plot label permutations
x = do.call('rbind', out)

for(i in unique(colnames(x))){
  y = as.data.frame(x[,i,drop=F])
  colnames(y) = 'v1'
  y$rank = rank(y$v1)
  y = y[order(y$v1,decreasing = T),]
  y$percentile = 1-(y$rank/nrow(y))
  fun = approxfun(x = y$v1, y = y$percentile)

  pval = fun(medians[[i]])

  plot = ggplot(y, aes(x = v1)) +
    geom_histogram(aes(y = ..density..),
                   colour = 1, fill = "white", bins = 250) +
   geom_vline(xintercept=medians[[i]], size=.5, color="black")+theme_classic() + annotate("text", x = -.2, y = 1, label = paste0('pval = ', pval)) + ggtitle(i)

  pdf(paste0('../plots/label_permutation_', i,'.pdf'), width = 5, height = 5)
  print(plot)
  dev.off()
}

print('done.')
saveRDS(out, '../data/other_analyses_outputs/permutation_analysis_label_perms.rds')
