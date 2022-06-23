library(ggplot2)
library(ggrepel)
cholest = readRDS('../data/other_analyses_outputs/cholesterol_analysis.rds')
cholest_genes = c(cholest$union_cholest_biosynth$gene, c("MYRF","MOG","PLP1","PLLP","MAG","OPALIN"))

# compare the pseudo bulk to the wilcox degs
oli_wilcox = readRDS('../data/differentially_expressed_genes_data/oli_wilcox_results.rds')
pseudo_bulk = read.csv('../data/other_analyses_outputs/pseudo_bulk_degs_single_cell_all_celltypes.csv')
nebula = readRDS('../data/other_analyses_outputs/nebula_oli_degs.rds')

merged = merge(pseudo_bulk[pseudo_bulk$celltype=='Oli',], oli_wilcox, by.y = 0, by.x = 'gene')
merged2 = merge(pseudo_bulk[pseudo_bulk$celltype=='Oli',], nebula, by = 'gene')
merged3 = merge(oli_wilcox, nebula, by.x = 0, by.y = 'gene')

pdf('../plots/comparison_wilcox_pseudobulk_degs_oli.pdf', width = 5, height = 5)
ggplot(merged, aes(x=logFC.x, y=auc)) +
geom_hex(bins = 100)+  geom_smooth(method=lm) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0.5)  + scale_fill_viridis_c()
dev.off()

pdf('../plots/comparison_pseudobulk_nebula.pdf', width = 7, height = 7)
ggplot(merged2, aes(x=logFC, y=logFC_APOE4, label = gene)) +
geom_hex(bins = 100) + geom_point(data = merged2[merged2$gene%in%cholest_genes,], col = 'red') + geom_smooth(method=lm) + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0)  + scale_fill_viridis_c() +
geom_text_repel(label = ifelse(merged2$gene%in%cholest_genes,merged2$gene,''), max.overlaps  = 100000)
dev.off()

# also show histogram of the effect sizes
pdf('../plots/pseudobulk_effect_sizes.pdf', width = 5, height = 5)
ggplot(merged2, aes(x = logFC)) +
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white", binwidth =0.01 ) +
  geom_density() + geom_vline(xintercept=(merged2[merged2$gene%in%c('DHCR24', 'FDFT1', 'OPALIN', 'PLLP'),'logFC']), size=.5, color="black")+theme_classic()
dev.off()

print(rank(merged2$logFC)[merged2$gene=='DHCR24'])/max(rank(merged2$logFC))
print(rank(merged2$logFC)[merged2$gene=='FDFT1'])/max(rank(merged2$logFC))
print(rank(merged2$logFC)[merged2$gene=='PLLP'])/max(rank(merged2$logFC))
print(rank(merged2$logFC)[merged2$gene=='OPALIN'])/max(rank(merged2$logFC))

pdf('../plots/comparison_wilcox_nebula.pdf', width = 5, height = 5)
ggplot(merged3, aes(x=auc, y=logFC_APOE4)) +
geom_hex(bins = 100)+  geom_smooth(method=lm) + theme_classic() + geom_vline(xintercept = 0.5) + geom_hline(yintercept = 0)  + scale_fill_viridis_c()
dev.off()
