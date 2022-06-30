########## APOE expression related to extended data figure 6 #############
###################################################################

# required libraries
library(ggplot2)

counts = read.csv('../data/iPSC_data/ipsc_bulk_no_drug.csv')

meta1 = read.table('../data/iPSC_data/ipsc_metadata.csv', sep = ',', header = F)
x = as.data.frame(read.table('../data/iPSC_data/ipsc_bulk_rnaseq_count_files/D19-13224-4367S_geneexp.txt', header = FALSE))
gene.df = as.data.frame(cbind(as.character(x$V1), as.character(x$V3)))
rownames(gene.df) = gene.df$V1
counts$gene = gene.df[counts[,1],'V2']
counts = counts[!duplicated(counts$gene),]
counts = counts[!is.na(counts$gene),]
rownames(counts) = counts$gene
counts$X = NULL
counts$gene = NULL

# get percentile in which APOE is expressed
means = rowMeans(counts)
df = as.data.frame(means)
df = df[!df$means==0,,drop = F]

print('APOE percentile expressed gene, among nonzero detected genes:')
d = df[order(df$means,decreasing = F),,drop = F]
percentile = (which(rownames(d) =='APOE')/nrow(d))*100

pdf('../plots/Extended_6/apoe_expression_ipsc.pdf', width = 5, height = 5)
ggplot(df, aes(x = (means))) +
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white", binwidth =0.5 ) +
  geom_density() + geom_vline(xintercept=(c((df['APOE',1]))), size=.5, color="black")+theme_classic() + xlim(c(0,50)) + annotate("text", x = 4, y=1, label = percentile)
dev.off()
