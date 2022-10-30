
library('readxl')
library(SingleCellExperiment)

## Supplementary Table S1.
# print('s1')
# sce = readRDS('../data/single_cell_data/single_cell_experiment_object_qced.rds')
# coldata = colData(sce)
# meta = coldata[!duplicated(coldata$projid),]
# names = c('projid','age_death', 'amyloid', 'braaksc', 'ceradsc', 'cogdx', 'msex', 'nft', 'pmi', 'apoe_genotype')
# write.csv(meta[,names], '../data/supplementary_tables/Supplementary_Table_S1.csv')

## Supplementary Table S2.
print('s2')
expressed = readRDS('../data/single_cell_data/expressed_genes_per_celltype.rds')
out = list()
for(i in names(expressed)){
  df = as.data.frame(expressed[[i]])
  df$celltype = i
  colnames(df) = c('expressed_gene', 'celltype')
  out[[i]] = df
}
out = do.call('rbind', out)
write.csv(out, '../data/supplementary_tables/Supplementary_Table_S2.csv')

## Supplementary Table S3.
print('s3')
av = readRDS('../data/single_cell_data/individual_level_averages_per_celltype.rds')
out = list()
for(i in names(av)){
  df = as.data.frame(av[[i]])
  df$celltype = i
  out[[i]] = df
}
out = do.call('rbind', out)
write.csv(out, '../data/supplementary_tables/Supplementary_Table_S3.csv')

## Supplementary Table S4
print('s4')
data = readRDS('../data/other_analyses_outputs/pathway_scores.rds')
go_bp = do.call('rbind', data$GO_BP$fits_all)

path_names_go = as.data.frame(read_excel('../data/supplementary_tables/all_path_names_go_renaming.xlsx', sheet = 'all_path_names_go_renaming'))
all_bp = merge(go_bp, path_names_go, by.x = 'names', by.y = 'name', all.x = T, all.y = T)
all_bp = all_bp[all_bp$P.Value<0.05,]
highlight = c('negative regulation of T cell receptor signaling pathway (GO:0050860)',
'positive regulation of response to cytokine stimulus (GO:0060760)',
'negative regulation of NIK/NF-kappaB signaling (GO:1901223)',
'positive regulation of tumor necrosis factor-mediated signaling pathway (GO:1903265)',
'I-kappaB kinase/NF-kappaB signaling (GO:0007249)',
'regulation of high voltage-gated calcium channel activity (GO:1901841)',
'regulation of voltage-gated calcium channel activity (GO:1901385)',
'positive regulation of tumor necrosis factor-mediated signaling pathway (GO:1903265)',
'positive regulation of excitatory postsynaptic potential (GO:2000463)',
'regulation of long-term neuronal synaptic plasticity (GO:0048169)',
'I-kappaB kinase/NF-kappaB signaling (GO:0007249)',
'ERK1 and ERK2 cascade (GO:0070371)',
'regulation of early endosome to late endosome transport (GO:2000641)',
'oligosaccharide-lipid intermediate biosynthetic process (GO:0006490)',
'negative regulation of amyloid-beta formation (GO:1902430)',
'DNA damage induced protein phosphorylation (GO:0006975)',
'intrinsic apoptotic signaling pathway in response to DNA damage by p53 class mediator (GO:0042771)',
'chaperone mediated protein folding requiring cofactor (GO:0051085)',
'chaperone-mediated protein complex assembly (GO:0051131)',
'regulation of cholesterol transport (GO:0032374)',
'cholesterol biosynthetic process (GO:0006695)',
'negative regulation of lipid storage (GO:0010888)',
'glycogen biosynthetic process (GO:0005978)',
'acetyl-CoA metabolic process (GO:0006084)')
all_bp$highlight = ifelse(all_bp$names%in%highlight, 'yes', 'no')
write.csv(all_bp, '../data/supplementary_tables/Supplementary_Table_S4.csv')

## Supplementary Table S5.
print('s5')
go_apoe = do.call('rbind', data$APOE$fits_all)
# add pathway renaming
path_names_apoe = as.data.frame(read_excel('../data/supplementary_tables/all_path_names_go_renaming.xlsx', sheet = 'apoe_terms'))
all_apoe = merge(go_apoe, path_names_apoe, by.x = 'names', by.y = 'full name', x.all = T, y.all = T)
write.csv(all_apoe, '../data/supplementary_tables/Supplementary_Table_S5.csv')

## Supplementary Table S6.
print('s6')
go_lipid = do.call('rbind', data$lipid$fits_all)
# add pathway renaming
path_names_lipids = as.data.frame(read_excel('../data/supplementary_tables/all_path_names_go_renaming.xlsx', sheet = 'lipid_terms'))
all_lipids = merge(go_lipid, path_names_lipids, by.x = 'names', by.y = 'full name', x.all = T, y.all = T)
write.csv(all_lipids, '../data/supplementary_tables/Supplementary_Table_S6.csv')

## Supplementary Table S7.
print('s7')
cc = read.csv('../data/supplementary_tables/cc_lipidomics_data.csv', check.names = F)
write.csv(cc, '../data/supplementary_tables/Supplementary_Table_S7.csv')

## Supplementary Table S8.
print('s8')
df = read.csv('../data/supplementary_tables/pfc_lipidomics_data.csv')
write.csv(df, '../data/supplementary_tables/Supplementary_Table_S8a.csv')
data_subset = read.csv('../data/supplementary_tables/pfc_lipidomics_data_metadata.csv')
write.csv(data_subset, '../data/supplementary_tables/Supplementary_Table_S8b.csv')
data_subset = read.csv('../data/supplementary_tables/pfc_lipidomics_qc_metrics.csv')
write.csv(data_subset, '../data/supplementary_tables/Supplementary_Table_S8c.csv')
data_subset = read.csv('../data/supplementary_tables/lipicomics_all_data_no_qc.csv')
write.csv(data_subset, '../data/supplementary_tables/Supplementary_Table_S8d.csv')

## Supplementary Table S9.
print('s9')
mat = read.csv('../data/supplementary_tables/ipsc_postmortem_merged_scaled_matrices_individual_level.csv')
write.csv(mat, '../data/supplementary_tables/Supplementary_Table_S9.csv')

## Supplementary Table S10.
print('s10')
norm_counts = read.csv('../data/iPSC_data/FPKM_table_OPC.txt', sep = '\t')
write.csv(norm_counts,'../data/supplementary_tables/Supplementary_Table_S10.csv')

## Supplementary Table S12.
print('s12')
degs = read.csv('../data/iPSC_data/OPC_DEG_statistics.txt', sep = '\t')
write.csv(degs, '../data/supplementary_tables/Supplementary_Table_S12.csv')

## Supplementary Table S13.
print('s13')
out_all = as.data.frame(readRDS('../data/differentially_expressed_genes_data/oli_wilcox_results.rds'))
out_all$grp = 'APOE34 & APOE44 vs APOE33 (AD and nonAD)'
oli_ad = as.data.frame(readRDS('../data/differentially_expressed_genes_data/oli_wilcox_results_AD.rds'))
oli_ad$grp = 'APOE34 vs APOE33 (AD only)'
oli_nonad = as.data.frame(readRDS('../data/differentially_expressed_genes_data/oli_wilcox_results_noAD.rds'))
oli_nonad$grp = 'APOE34 vs APOE33 (nonAD only)'
all_data = rbind(out_all, oli_ad, oli_nonad)
write.csv(all_data, '../data/supplementary_tables/Supplementary_Table_S13.csv')

## Supplementary Table S14.
print('s14')
nebula = readRDS('../data/differentially_expressed_genes_data/E4_nebula_associations_by_celltype.rds')$Oli
nebula$padj = p.adjust(nebula$p_Apoe_e4yes, 'fdr')
write.csv(nebula, '../data/supplementary_tables/Supplementary_Table_S14.csv')

## Supplementary Table S15.
print('s15')
degs = read.csv('../data/other_analyses_outputs/pseudo_bulk_degs_single_cell_all_celltypes.csv')
write.csv(degs, '../data/supplementary_tables/Supplementary_Table_S15.csv')
