########## fgsea analysis for extended data figure 2 ##########
###############################################################

# required packages
library(fgsea)
print('running fgsea analysis for extended data figure 2')

neb = readRDS('../data/differentially_expressed_genes_data/E4_nebula_associations_by_celltype.rds')

expressed = readRDS('../data/single_cell_data/expressed_genes_per_celltype.rds')
scores = sign(neb$Oli$logFC_Apoe_e4yes) * -log10(neb$Oli$p_Apoe_e4yes)
names(scores) = rownames(neb$Oli)
sorted = sort(scores)

pathways = readRDS('../data/other_analyses_outputs/pathways.rds')

low_removed_lipid_paths = pathways$pathways$low_removed_lipid_associated

fgseaRes <- fgsea(pathways = low_removed_lipid_paths$Oli, stats = sorted,nperm = 10000, minSize = 5, maxSize = 500)

df = fgseaRes[fgseaRes$pval<0.05,]
x = df$NES
names(x) = df$pathway

pdf('../plots/Extended_2/fgsea_lipid_associated_paths.pdf', width = 4, height = 4)
barplot(x[order(x)], las = 1, horiz = T,xlab = 'normalized enrichment score')
dev.off()

print('done.')
