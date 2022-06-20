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
# TODO: clean up the data files for upload
print('loading the data..')
e = readRDS('../../submission_code_09012021/data/Summary.DE.celltype.rds')
expressed = metadata(e)$expressed.genes
av_expression = readRDS('../../submission_code_09012021/data/Averages.by.celltype.by.individual.rds')
summary = as.data.frame(read_excel('../../submission_code_09012021/data/Overview_PFC_projids_data_384.xlsx'))
summary = summary[!duplicated(summary[,'projid...2']),]
rownames(summary) = summary[['projid...2']]

# load all the pathway terms
print('loading the pathway terms...')
bp = GSA.read.gmt('../data/GO_Biological_Process_2018.txt')
bp_g = bp$genesets
names(bp_g) = bp$geneset.names

cyc = GSA.read.gmt('../data/HumanCyc_2016.txt')
cyc_g = cyc$genesets
names(cyc_g) = cyc$geneset.names

kegg = GSA.read.gmt('../data/KEGG_2019_Human.txt')
kegg_g = kegg$genesets
names(kegg_g) = kegg$geneset.names

reactome = GSA.read.gmt('../data/Reactome_2016.txt')
reactome_g = reactome$genesets
names(reactome_g) = reactome$geneset.names

all_paths = c(cyc_g, reactome_g, kegg_g, bp_g)

# filter out the lowly expressed genes (10%, for each celltype)
print('filtering pathways by expression...')
low_removed_bp = filter_lowly_exp_genes(expressed, bp_g)
low_removed_reactome = filter_lowly_exp_genes(expressed, reactome_g)
low_removed_kegg = filter_lowly_exp_genes(expressed, kegg_g)
low_removed_humancyc = filter_lowly_exp_genes(expressed, cyc_g)
low_removed_all_terms = filter_lowly_exp_genes(expressed, all_paths)

# add the filtered pathways to the dataset to be saved
all_data[['pathways']][['low_removed_reactome']] = low_removed_reactome
all_data[['pathways']][['low_removed_kegg']] = low_removed_kegg
all_data[['pathways']][['low_removed_humancyc']] = low_removed_humancyc
all_data[['pathways']][['low_removed_bp']] = low_removed_bp

all_data[['pathways']][['all_reactome']] = reactome_g
all_data[['pathways']][['all_kegg']] = kegg_g
all_data[['pathways']][['all_humancyc']] = cyc_g
all_data[['pathways']][['all_bp']] = bp_g
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
saveRDS(all_data, '../data_outputs/pathways.rds')
