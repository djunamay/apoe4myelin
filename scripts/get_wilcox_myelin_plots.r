
# show the degs for myelin genes with and without AD
ad = readRDS('../data/differentially_expressed_genes_data/oli_wilcox_results_AD.rds')
myelination <- c("MYRF","MOG","PLP1","PLLP","MAG","OPALIN")
gg = ad[myelination,]
gg = gg[gg$padj<0.05,]
ggg = gg$logFC
names(ggg) = rownames(gg)

pdf('../plots/wilcox_ad_myelin.pdf', width = 3, height = 3.5)
barplot(ggg[order(ggg, decreasing = T)], las = 2, horiz = T, main = 'AD')
dev.off()

# show the degs for myelin genes with and without AD
no = readRDS('../data/differentially_expressed_genes_data/oli_wilcox_results_noAD.rds')
myelination <- c("MYRF","MOG","PLP1","PLLP","MAG","OPALIN")
gg = no[myelination,]
gg = gg[gg$padj<0.05,]
ggg = gg$logFC
names(ggg) = rownames(gg)

pdf('../plots/wilcox_noad_myelin.pdf', width = 3, height = 3.5)
barplot(ggg[order(ggg, decreasing = T)], las = 2, horiz = T, main = 'no AD')
dev.off()

print('done.')
