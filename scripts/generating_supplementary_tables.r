
library('readxl')

## differential pathway results BP
data = readRDS('../data/other_analyses_outputs/pathway_scores.rds')
go_bp = do.call('rbind', data$GO_BP$fits_all)
# add pathway renaming
path_names_go = as.data.frame(read_excel('../data/supplementary_tables/all_path_names_go_renaming.xlsx', sheet = 'all_path_names_go_renaming'))
all_bp = merge(go_bp, path_names_go, by.x = 'names', by.y = 'name', keep_all = T)
write.csv(all_bp, '../data/supplementary_tables/Supplementary_Table_S4.csv')

## differential pathway results APOE
go_apoe = do.call('rbind', data$APOE$fits_all)
# add pathway renaming
path_names_apoe = as.data.frame(read_excel('../data/supplementary_tables/all_path_names_go_renaming.xlsx', sheet = 'apoe_terms'))
all_apoe = merge(go_apoe, path_names_apoe, by.x = 'names', by.y = 'full name', keep_all = T)
write.csv(all_apoe, '../data/supplementary_tables/Supplementary_Table_S5.csv')

## differential pathway results lipids
go_lipid = do.call('rbind', data$lipid$fits_all)
# add pathway renaming
path_names_lipids = as.data.frame(read_excel('../data/supplementary_tables/all_path_names_go_renaming.xlsx', sheet = 'lipid_terms'))
all_lipids = merge(go_lipid, path_names_lipids, by.x = 'names', by.y = 'full name', keep_all = T)
write.csv(all_lipids, '../data/supplementary_tables/Supplementary_Table_S6.csv')

## CC lipidomics data
cc = read.csv('../data/supplementary_tables/cc_lipidomics_data.csv', check.names = F)
write.csv(cc, '../data/supplementary_tables/Supplementary_Table_S7.csv')

## PFC lipidomics data

## ipsc degs
degs = read.csv('../data/iPSC_data/apoe3_v_apoe4_degs.csv')
write.csv(degs, '../data/supplementary_tables/Supplementary_Table_S12.csv')

## ipsc counts
out_APOE4 = read.table('../data/iPSC_data/apoe3_v_apoe4_degs.csv', sep = ',', header = T)
norm_counts = read.csv('../data/iPSC_data/ipsc_bulk_no_drug.csv')
rownames(norm_counts) = norm_counts$X
norm_counts$X = NULL
rownames(out_APOE4) = out_APOE4$X
out_APOE4 = out_APOE4[!duplicated(out_APOE4$res.gene_name),]
shared = intersect(rownames(out_APOE4), rownames(norm_counts))
norm_counts = norm_counts[shared,]
rownames(norm_counts) = out_APOE4[shared,'res.gene_name']

meta1 = read.table('../data/iPSC_data/ipsc_metadata.csv', sep = ',', header = F)
var = as.character(meta1$V2)
colnames(norm_counts) = var
write.csv(norm_counts,'../data/supplementary_tables/Supplementary_Table_S10.csv')

## pseudo bulk degs
degs = read.csv('../data/other_analyses_outputs/pseudo_bulk_degs_single_cell_all_celltypes.csv')
write.csv(degs, '../data/supplementary_tables/Supplementary_Table_S15.csv')

## wilcox degs for oligodendrocytes
out_all = as.data.frame(readRDS('../data/differentially_expressed_genes_data/oli_wilcox_results.rds'))
out_all$grp = 'APOE34 & APOE44 vs APOE33 (AD and nonAD)'
oli_ad = as.data.frame(readRDS('../data/differentially_expressed_genes_data/oli_wilcox_results_AD.rds'))
oli_ad$grp = 'APOE34 vs APOE33 (AD only)'
oli_nonad = as.data.frame(readRDS('../data/differentially_expressed_genes_data/oli_wilcox_results_noAD.rds'))
oli_nonad$grp = 'APOE34 vs APOE33 (nonAD only)'
all_data = rbind(out_all, oli_ad, oli_nonad)
write.csv(all_data, '../data/supplementary_tables/Supplementary_Table_S13.csv')

## nebula degs for oligodendrocytes
data = readRDS('../data/other_analyses_outputs/nebula_oli_degs.rds')
write.csv(data, '../data/supplementary_tables/Supplementary_Table_S14.csv')
