# get the data
counts = read.csv('../data_outputs/ipsc_bulk_no_drug.csv')

meta1 = read.table('../raw_data/APOE_ipsc_opc_sequencing_meta.csv', sep = ',', header = F)
x = as.data.frame(read.table('../raw_data/191210Tsa/RNAseqQC-4367S/D19-13224-4367S_geneexp.txt', header = FALSE))
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
print(which(rownames(d) =='APOE')/nrow(d))

print('fraction of genes captured in the untransformed histogram:')
x = table(df$means<50)
print(x[2]/sum(x))

library(ggplot2)

pdf('../plots/apoe_expression_ipsc.pdf', width = 5, height = 5)
ggplot(df, aes(x = (means))) +
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white", binwidth =0.5 ) +
  geom_density() + geom_vline(xintercept=(c((df['APOE',1]))), size=.5, color="black")+theme_classic() + xlim(c(0,50))
dev.off()
