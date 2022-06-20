########## script 1 in APOE4 impairs myelination via cholesterol dysregulation in oligodendrocytes ##########
#############################################################################################################
#TODO: remove the genes that are empty from all the genesets

##### required libraries ####
source('../functions/pathway_analyses.r')
library('readxl')
library('SingleCellExperiment')
library('GSA')

all_data = list()

# load the data
print('loading the data..')
expressed = readRDS('../data/single_cell_data/expressed_genes_per_celltype.rds')
av_expression = readRDS('../data/single_cell_data/individual_level_averages_per_celltype.rds')
summary = read.csv('../data/single_cell_data/metadata_by_individual.csv')

# load all the pathway terms
print('loading the pathway terms...')
bp = read.geneset('../data/pathway_databases/GO_Biological_Process_2018.txt')
cy = read.geneset('../data/pathway_databases/HumanCyc_2016.txt')
ke = read.geneset('../data/pathway_databases/KEGG_2019_Human.txt')
re = read.geneset('../data/pathway_databases/Reactome_2016.txt')

all_paths = c(bp, cy, ke, re)

# filter out the lowly expressed genes (10%, for each celltype)
print('filtering pathways by expression...')
low_removed_bp = filter_lowly_exp_genes(expressed, bp)
low_removed_reactome = filter_lowly_exp_genes(expressed, re)
low_removed_kegg = filter_lowly_exp_genes(expressed, ke)
low_removed_humancyc = filter_lowly_exp_genes(expressed, cy)
low_removed_all_terms = filter_lowly_exp_genes(expressed, all_paths)

# add the filtered pathways to the dataset to be saved
all_data[['pathways']][['low_removed_reactome']] = low_removed_reactome
all_data[['pathways']][['low_removed_kegg']] = low_removed_kegg
all_data[['pathways']][['low_removed_humancyc']] = low_removed_humancyc
all_data[['pathways']][['low_removed_bp']] = low_removed_bp

all_data[['pathways']][['all_reactome']] = re
all_data[['pathways']][['all_kegg']] = ke
all_data[['pathways']][['all_humancyc']] = cy
all_data[['pathways']][['all_bp']] = bp
all_data[['pathways']][['all']] = all_paths

# get apoe-associated genesets
print('filtering APOE4 terms...')
gene = 'APOE'
index = unlist(lapply(names(all_paths), function(x) c(gene)%in%all_paths[[x]]))
apoe_gsets = all_paths[index]
gs = lapply(names(apoe_gsets), function(x) apoe_gsets[[x]][apoe_gsets[[x]]!=''])
names(gs) = names(apoe_gsets)
apoe_gsets = gs
apoe_gsets_low_removed = filter_lowly_exp_genes(expressed, apoe_gsets)
all_data[['pathways']][['apoe_gsets_low_removed']] = apoe_gsets_low_removed
all_data[['pathways']][['apoe_gsets_all']] = apoe_gsets

# get lipid-associated genesets
lipid_keys = c('sterol','athero','cholest','LDL','HDL','lipoprotein','triglyceride','TAG','DAG','lipid','steroid','fatty acid', 'ceramide')
lipid_paths = get_gset_names_by_category(lipid_keys, names(all_paths))
low_removed_lipid_paths = filter_lowly_exp_genes(expressed, all_paths[lipid_paths])
all_data[['pathways']][['all_lipid_associated']] = lipid_paths
all_data[['pathways']][['low_removed_lipid_associated']] = low_removed_lipid_paths

# save all the data
saveRDS(all_data, '../data/other_analyses_outputs/pathways.rds')
print('done.')
