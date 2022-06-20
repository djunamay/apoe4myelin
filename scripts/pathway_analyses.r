########## script 2 in APOE4 impairs myelination via cholesterol dysregulation in oligodendrocytes ##########
#############################################################################################################

##### required libraries ####
source('../functions/pathway_analyses.r')
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
pathways = readRDS('../data_outputs/pathways.rds')
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

##############################
###### GO BP PATHWAYS #####
##############################
print('analysis on GO BP pathways...')
# run GSVA on GO pathways
print('running GSVA...')
order = c('In', 'Ex', 'Oli', 'Opc', 'Mic','Ast')
out_bp_terms = lapply(order, function(x) t(gsva(as.matrix(av_expression[[x]]), low_removed_bp[[x]], mx.diff=TRUE, verbose=F, kcdf=c("Gaussian"), min.sz=5, max.sz = 150, parallel.sz=16)))
names(out_bp_terms) = order
all_data[['res']][['BP']][['gsva_out']] = out_bp_terms

# get linear model fits
print('getting linear model fits...')
summary$APOE4 = ifelse(summary$apoe_genotype == '33','e3','e4')
fits = get_fits(out_bp_terms, summary)
all_data[['res']][['BP']][['fits']] = fits

# get matrix of scores for heatmap
print('visualizing scores as heatmap...')
scores = get_scores(fits)
all_data[['res']][['GO_BP']][['scores']] = scores
names = unique(unname(unlist(lapply(names(scores$all), function(x) rownames(scores$all[[x]])))))
mat = get_matrix(scores$all, names)
mat = mat[unname(rowSums(abs(mat)>1.3)>0),]
index = unique(unname(unlist(lapply(colnames(mat), function(x) order(abs(mat[[x]]),decreasing = T)[1:20]))))
mat = mat[index,]
all_data[['res']][['GO_BP']][['scores_top_20']] = mat

#####################################
###### APOE-ASSOCIATED PATHWAYS #####
#####################################
print('analysis on APOE-associated pathways...')
# run GSVA on APOE-associated pathways
print('running GSVA...')
order = c('In', 'Ex', 'Oli', 'Opc', 'Mic','Ast')
out_apoe_terms = lapply(order, function(x) t(gsva(as.matrix(av_expression[[x]]), apoe_gsets_low_removed[[x]], mx.diff=TRUE, verbose=F, kcdf=c("Gaussian"), min.sz=5, max.sz = 150, parallel.sz=16)))
names(out_apoe_terms) = order
all_data[['res']][['apoe_associated']][['gsva_out']] = out_apoe_terms

# get linear model fits
print('getting linear model fits...')
summary$APOE4 = ifelse(summary$apoe_genotype == '33','e3','e4')
fits = get_fits(out_apoe_terms, summary) # TODO make the trace more informative
all_data[['res']][['apoe_associated']][['fits']] = fits

# visualize scores as heatmap
print('getting scores...')
scores = get_scores(fits)
all_data[['res']][['apoe_associated']][['scores']] = scores

# show the heatmap for the APOE-associated pathways
print('drawing the apoe-associated heatmap...')
names = unique(unname(unlist(lapply(names(scores$all), function(x) rownames(scores$all[[x]])))))
mat = get_matrix(scores$all, names)
df = mat[unname(rowSums(abs(mat)>1.3)>0),]
all_data[['res']][['apoe_associated']][['scores_sig']] = df

#####################################
###### LIPID-ASSOCIATED PATHWAYS #####
#####################################
print('analysis on lipid-associated pathways...')
# run GSVA on reactome pathways
print('running GSVA...')
order = c('In', 'Ex', 'Oli', 'Opc', 'Mic','Ast')
out_lipid_terms = lapply(order, function(x) t(gsva(as.matrix(av_expression[[x]]), low_removed_lipid_paths[[x]], mx.diff=TRUE, verbose=F, kcdf=c("Gaussian"), min.sz=5, max.sz = 150, parallel.sz=16)))
names(out_lipid_terms) = order
all_data[['res']][['lipid_associated']][['gsva_out']] = out_lipid_terms

# get linear model fits
print('getting linear model fits...')
summary$APOE4 = ifelse(summary$apoe_genotype == '33','e3','e4')
fits = get_fits(out_lipid_terms, summary)
all_data[['res']][['lipid_associated']][['fits']] = fits

# visualize scores as heatmap
print('getting scores...')
scores = get_scores(fits)
all_data[['res']][['lipid_associated']][['scores']] = scores

# show the heatmap
print('drawing the apoe-associated heatmap...')
names = unique(unname(unlist(lapply(names(scores$all), function(x) rownames(scores$all[[x]])))))
mat = get_matrix(scores$all, names)
df = mat[unname(rowSums(abs(mat)>1.3)>0),]
all_data[['res']][['lipid_associated']][['scores_sig']] = df

# save all the data
saveRDS(all_data, '../data_outputs/pathway_scores.rds')
print('done.')
