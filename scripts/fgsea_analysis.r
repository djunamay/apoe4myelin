# adding on June 9
# additional analysis of nebula rankings to test for cholesterol dysregulation
# add to methods in paper all the details about pathway estimation
neb = readRDS('../data/nebula.E4.associations.rds')

library(ACTIONet)
go = readRDS('/home/djuna/Documents/data/GO.rds')

scores = neb$nebula.out.processed.APOE4.scores
enrich = ACTIONet::assess.geneset.enrichment.from.scores(scores, go)
head(enrich$logPvals)
head(enrich$scores)

colnames(enrich$logPvals) = colnames(scores)

pvals = as.data.frame(exp(-1*(enrich$logPvals)))
pvals$padj.oli = unname(p.adjust(pvals[,'Oli_Apoe_e4yes'],'fdr'))

pvals = pvals[order(pvals$padj.oli,decreasing=F),]
head(pvals)

e = readRDS('../data/data/Summary.DE.celltype.rds')
expressed = metadata(e)$expressed.genes
sorted = sort(scores[,'Oli_Apoe_e4yes'])
sorted = sorted[names(sorted)%in%expressed$Oli]

fgseaRes <- fgsea(pathways = go, stats = sorted,nperm = 10000, minSize = 15, maxSize = 500)
fgseaRes[order(fgseaRes$pval,decreasing = F),]


df = (neb$GSEA.out$Oli_Apoe_e4yes)
df = df[df$padj<0.05,]
paths = c('cholesterol biosynthetic process (GO:0006695)',
          'cholesterol metabolic process (GO:0008203)',
          'cholesterol homeostasis (GO:0042632)',
          'sterol biosynthetic process (GO:0016126)',
          'sterol homeostasis (GO:0055092)'


)
f = df[df$pathway%in%paths,c('pathway','NES')]
ff = f$NES
names(ff) = f$pathway

#pdf('../results/nbmm_oli_cholest.pdf',width = 2, height = 3)
barplot(ff, las = 1, horiz = T, xlab = 'normalized enrichment score')
#dev.off()
