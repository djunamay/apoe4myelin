# load the wilcox degs
degs = readRDS('../data/oli_wilcox_results.rds')

options(repr.plot.width = 6, repr.plot.height =3.5, repr.plot.res = 100)
# these are a selection of genes chosen by Joel for figure 2
genes = c('PLLP', 'MYRF', 'MAG', 'OPALIN', 'SREBF1', 'PLP1', 'SCAP',
         'MVK', 'FDPS', 'ABCG1', 'HMGCS1', 'IDI1', 'LDLR', 'INSIG1', 'SREBF2', 'SQLE', 'DHCR7','DHCR24', 'FDFT1', 'LSS')
gg = degs[genes,]
gg = gg[gg$padj<0.05,]
ggg = gg$logFC
names(ggg) = rownames(gg)
#pdf('../figure_outputs/wilcox_all.pdf', width = 8, height = 3.5)
barplot(ggg, las = 2)
#dev.off()

# show the degs for myelin genes with and without AD
options(repr.plot.width = 3, repr.plot.height =3.5, repr.plot.res = 100)

ad = readRDS('../data/oli_wilcox_results_AD.rds')
myelination <- c("MYRF","MOG","PLP1","PLLP","MAG","OPALIN")
gg = ad[myelination,]
gg = gg[gg$padj<0.05,]
ggg = gg$logFC
names(ggg) = rownames(gg)

#pdf('../figure_outputs/wilcox_ad.pdf', width = 3, height = 3.5)
barplot(ggg[order(ggg, decreasing = T)], las = 2, horiz = T, main = 'AD')
#dev.off()

# show the degs for myelin genes with and without AD
options(repr.plot.width = 3, repr.plot.height =3.5, repr.plot.res = 100)

no = readRDS('../data/oli_wilcox_results_noAD.rds')
myelination <- c("MYRF","MOG","PLP1","PLLP","MAG","OPALIN")
gg = no[myelination,]
gg = gg[gg$padj<0.05,]
ggg = gg$logFC
names(ggg) = rownames(gg)

#pdf('../figure_outputs/wilcox_noad.pdf', width = 3, height = 3.5)
barplot(ggg[order(ggg, decreasing = T)], las = 2, horiz = T, main = 'no AD')
#dev.off()
