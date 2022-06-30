########## iPSC analysis related to extended data figure 8 #############
##########################################################################

# required packages
library(GSVA)
library(limma)
library(tidyr)
library(reshape2)
library(ggplot2)
source('../functions/pathway_analyses.r')

out_APOE4 = read.table('../data/iPSC_data/OPC_DEG_statistics.txt',  header = T)
pathways = readRDS('../data/other_analyses_outputs/pathways.rds')

# look at the cholesterol biosynthesis signature
data = readRDS('../data/other_analyses_outputs/cholesterol_analysis.rds')
biosynth_genes = data[['union_cholest_biosynth']]$genes

out = list()
out[['cholest_bio']] = biosynth_genes
out[['ATF6']] = c('HSP90B1', 'CALR', 'HSPA5', 'MBTPS1', 'ATF6')

curr = out_APOE4
curr = curr[!duplicated(curr$gene_id),]
rownames(curr) = curr$gene_id

l = list()
for(i in names(out)){
    df = na.omit(curr[out[[i]],])
    df$score = df$log2.fold_change.#) * -log10(df$q_value)
    df$grp = i
    l[[i]] = df[df$q_value<0.05,]
}

for(f in names(l)){
    x = l[[f]]$score
    names(x) = rownames(l[[f]])
    pdf(paste0('../plots/Extended_8/', f,'.pdf'), width = 3, height = 5)
    print(barplot(x[order(x)], las = 1, horiz = T))
    dev.off()
}

# also show myelin genes
genes = c('PLP1','OPALIN','PLLP','MYRF','MAG','MOG')
df = curr[genes,]
df$score = sign(df$log2.fold_change.) * -log10(df$q_value)
x = df$score
names(x) = rownames(df)

pdf('../plots/myelin_ipsc.pdf', width = 3, height = 3)
print(barplot(x[order(x)], las = 1, horiz = T))
dev.off()

# show lipid/cholesterol-related pathway enrichment
all_paths = pathways$pathways$all
norm_counts =  read.table('../data/iPSC_data/FPKM_table_OPC.txt', header = TRUE)
rownames(norm_counts) = norm_counts$gene
norm_counts$gene = NULL
rownames(out_APOE4) = out_APOE4$gene
shared = intersect(rownames(out_APOE4), rownames(norm_counts))
norm_counts = norm_counts[shared,]
rownames(norm_counts) = out_APOE4[shared,'gene']

# get the metadata
var = as.data.frame(ifelse(startsWith(colnames(norm_counts), 'E3'), 0, 1))
colnames(var) = c('APOE')

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
pdf('../plots/Extended_8/upr_paths_ips.pdf', width = 3, height = 3)
print(barplot(x[order(x)], las = 1, horiz = T))
dev.off()

x = allgenesets_cholesterol[allgenesets_cholesterol$P.Value<0.05,'logFC']
names(x) = rownames(allgenesets_cholesterol[allgenesets_cholesterol$P.Value<0.05,])
pdf('../plots/Extended_8/cholest_paths_ips.pdf', width = 3, height = 3)
print(barplot(x[order(x)], las = 1, horiz = T))
dev.off()

# show other key genes of interest
colnames(norm_counts) = as.character(meta1$V2)
colnames(norm_counts) = paste0(colnames(norm_counts),'_', c(1,2,3))
df = na.omit(norm_counts[c('SOAT1', 'SOAT2', 'CYP46A1'),])
df$gene = rownames(df)
x = melt(df)
x$variable = ifelse(startsWith(as.character(x$variable), 'APOE3'), 'APOE3', 'APOE4')
pdf('../plots/Extended_8/SOAT1_CYP_boxplots_ipsc.pdf', width = 3, height = 4)
ggplot(x, aes(x=variable, y=value, col = variable)) +
      geom_boxplot(width = .5) + geom_jitter() + facet_wrap(. ~ gene,  scales="free_y", nrow = 2) + theme_classic()
dev.off()
print('done.')


out_APOE4[out_APOE4$gene%in%c('SOAT1', 'CYP46A1'),]
