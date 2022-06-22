########## extended script 2 in APOE4 impairs myelination via cholesterol dysregulation in oligodendrocytes ##########
#############################################################################################################

library(GSVA)
library(limma)
library(tidyr)
source('../functions/pathway_analyses.r')

out_APOE4 = read.table('../data_outputs/apoe3_v_apoe4_degs.csv', sep = ',', header = T)
pathways = readRDS('../data_outputs/pathways.rds')

# look at the cholesterol biosynthesis signature
f1_data = readRDS('../data_outputs/pathway_scores.rds')
oli.lipid.fits = f1_data$res$lipid_associated$fits$Oli
lipid.paths.oli = pathways$pathways$low_removed_lipid_associated$Oli
paths = rownames(oli.lipid.fits[oli.lipid.fits$P.Value<0.05,])
biosynth_genes = unique(unname(unlist(lipid.paths.oli[paths])))

out = list()
out[['cholest_bio']] = biosynth_genes
out[['ATF6']] = c('HSP90B1', 'CALR', 'HSPA5', 'MBTPS1', 'ATF6')

curr = out_APOE4
curr = curr[!duplicated(curr$res.gene_name),]
rownames(curr) = curr$res.gene_name

l = list()
for(i in names(out)){
    df = na.omit(curr[out[[i]],])
    df$score = sign(df$res.logFC) * -log10(df$res.adj.P.Val)
    df$grp = i
    l[[i]] = df[df$res.adj.P.Val<0.05,]
}

for(f in names(l)){
    x = l[[f]]$score
    names(x) = rownames(l[[f]])
    pdf(paste0('../plots/', f,'.pdf'), width = 3, height = 3)
    print(barplot(x[order(x)], las = 1, horiz = T))
    dev.off()
}

# also show myelin genes
genes = c('PLP1','OPALIN','PLLP','MYRF','MAG','MOG')
df = curr[genes,]
df$score = sign(df$res.logFC) * -log10(df$res.adj.P.Val)
x = df$score
names(x) = rownames(df)

pdf('../plots/myelin_ipsc.pdf', width = 3, height = 3)
print(barplot(x[order(x)], las = 1, horiz = T))
dev.off()

# show lipid/cholesterol-related pathway enrichment
all_paths = pathways$pathways$all
norm_counts =  read.csv('../data_outputs/ipsc_bulk_no_drug.csv')
rownames(norm_counts) = norm_counts$X
norm_counts$X = NULL
rownames(out_APOE4) = out_APOE4$X
out_APOE4 = out_APOE4[!duplicated(out_APOE4$res.gene_name),]
shared = intersect(rownames(out_APOE4), rownames(norm_counts))
norm_counts = norm_counts[shared,]
rownames(norm_counts) = out_APOE4[shared,'res.gene_name']

# get the metadata
meta1 = read.table('../raw_data/APOE_ipsc_opc_sequencing_meta.csv', sep = ',', header = F)
var = as.character(meta1$V2)
var = as.data.frame(ifelse(var=='APOE4', 1, 0))
colnames(var) = 'APOE'

# get cholesterol-associated pathways
lipid_keys = c('sterol','steroid')
lipid_paths = get_gset_names_by_category(lipid_keys, names(all_paths))
out = gsva(as.matrix(norm_counts), all_paths[lipid_paths], mx.diff=TRUE, verbose=F, kcdf=c("Gaussian"), min.sz=5, max.sz = 150, parallel.sz=16)
mod = model.matrix(~APOE, data=var)
fit = lmFit(out, design=mod)
fit = eBayes(fit)
allgenesets_cholesterol = topTable(fit, coef='APOE', number=Inf, confint = T) %>% .[order(.$P.Value, decreasing = F),]

# get ATF/UPR-associated pathways
upr_keys = c('ATF', 'unfolded protein')
upr_paths = get_gset_names_by_category(upr_keys, names(all_paths))
out = gsva(as.matrix(norm_counts), all_paths[upr_paths], mx.diff=TRUE, verbose=F, kcdf=c("Gaussian"), min.sz=5, max.sz = 150, parallel.sz=16)
mod = model.matrix(~APOE, data=var)
fit = lmFit(out, design=mod)
fit = eBayes(fit)
allgenesets_upr = topTable(fit, coef='APOE', number=Inf, confint = T) %>% .[order(.$P.Value, decreasing = F),]

# show the pathways as barplots
x = allgenesets_upr[allgenesets_upr$P.Value<0.05,'logFC']
names(x) = rownames(allgenesets_upr[allgenesets_upr$P.Value<0.05,])
pdf('../plots/upr_paths_ips.pdf', width = 3, height = 3)
print(barplot(x[order(x)], las = 1, horiz = T))
dev.off()

x = allgenesets_cholesterol[allgenesets_cholesterol$P.Value<0.05,'logFC']
names(x) = rownames(allgenesets_cholesterol[allgenesets_cholesterol$P.Value<0.05,])
pdf('../plots/cholest_paths_ips.pdf', width = 3, height = 3)
print(barplot(x[order(x)], las = 1, horiz = T))
dev.off()
