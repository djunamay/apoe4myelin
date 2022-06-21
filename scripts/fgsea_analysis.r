# adding on June 9
# additional analysis of nebula rankings to test for cholesterol dysregulation
# add to methods in paper all the details about pathway estimation
library(fgsea)
neb = readRDS('../data/differentially_expressed_genes_data/E4_nebula_associations_by_celltype.rds')
set.seed(5)
# library(ACTIONet)
# go = readRDS('/home/djuna/Documents/data/GO.rds')
#
# scores = neb$nebula.out.processed.APOE4.scores
# enrich = ACTIONet::assess.geneset.enrichment.from.scores(scores, go)
# head(enrich$logPvals)
# head(enrich$scores)
#
# colnames(enrich$logPvals) = colnames(scores)
#
# pvals = as.data.frame(exp(-1*(enrich$logPvals)))
# pvals$padj.oli = unname(p.adjust(pvals[,'Oli_Apoe_e4yes'],'fdr'))
#
# pvals = pvals[order(pvals$padj.oli,decreasing=F),]
# head(pvals)

expressed = readRDS('../data/single_cell_data/expressed_genes_per_celltype.rds')
scores = sign(neb$Oli$logFC_Apoe_e4yes) * -log10(neb$Oli$p_Apoe_e4yes)
names(scores) = rownames(neb$Oli)
sorted = sort(scores)
#sorted = sorted[names(sorted)%in%expressed$Oli]

pathways = readRDS('../data/other_analyses_outputs/pathways.rds')

low_removed_lipid_paths = pathways$pathways$low_removed_lipid_associated

fgseaRes <- fgsea(pathways = low_removed_lipid_paths$Oli, stats = sorted,nperm = 10000, minSize = 5, maxSize = 500)

df = fgseaRes[fgseaRes$pval<0.05,]
x = df$NES
names(x) = df$pathway
pdf('../plots/fgsea_lipid_associated_paths.pdf', width = 4, height = 4)
barplot(x[order(x)], las = 1, horiz = T,xlab = 'normalized enrichment score')
dev.off()
