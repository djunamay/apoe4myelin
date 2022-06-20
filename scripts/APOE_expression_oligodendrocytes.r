# population of APOE-expressing oligodendrocytes

# load sce
sce = readRDS('../../submission_code_09012021/data/snRNA.data.sce.rds')
oli = sce[,sce$cell.type=='Oli']

# get APOE expression
d = as.data.frame(logcounts(oli)['APOE',])
colnames(d) = c('APOE')

# how is APOE detected relative to other genes?
library(ggplot2)

cts = counts(oli)
cts_binary = cts>0
sums = rowSums(cts_binary)
det.rate = sums/ncol(cts_binary)
df = as.data.frame(det.rate)

options(repr.plot.width = 3, repr.plot.height = 3, repr.plot.res = 300)

print('order of genes=')
print(c('APOE', 'NEUROD1', 'OLIG2', 'MBP'))
print(c(df['APOE',1], df['NEUROD1',1], df['OLIG2',1],df['MBP',1]))

pdf('../plots/distribution_apoe.pdf', width = 3, height = 3)
ggplot(df, aes(x = det.rate)) +
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white", binwidth =0.01 ) +
  geom_density() + geom_vline(xintercept=(c(df['APOE',1], df['NEUROD1',1], df['OLIG2',1],df['MBP',1])), size=.5, color="black")+theme_classic()
dev.off()

# how large is the population of APOE-detecting Oligodendrocytes

d = as.data.frame(logcounts(oli)['APOE',])
colnames(d) = c('APOE')
clust = kmeans(d$APOE, 2)
clust = (clust$cluster)

d$clust = clust
pdf('../plots/Apoe_density.pdf', width = 3, height = 3)
ggplot(d, aes(x = APOE)) +
  geom_histogram(
                 colour = 1, fill = "white", binwidth =0.1 ) +
  geom_density(aes(fill = clust, alpha = 0.5)) +theme_classic() +facet_wrap( ~ clust, nrow = 3, scales = 'free_y')
dev.off()
