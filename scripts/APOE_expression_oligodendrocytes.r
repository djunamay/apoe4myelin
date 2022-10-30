########## APOE expression related to extended data figure 6 #############
###################################################################

print('|| plotting APOE expression in oligodendrocytes (post-mortem)... ||')

# required packages
library(ggplot2)
library(SingleCellExperiment)
# load sce
print('loading')
sce = readRDS('../data/single_cell_data/single_cell_experiment_object_qced.rds')
print('subsetting')
oli = sce[,sce$cell.type=='Oli']

# get APOE expression
print('apoe expression')
d = as.data.frame(logcounts(oli)['APOE',])
colnames(d) = c('APOE')

# how is APOE detected relative to other genes?
print('apoe detection')
cts = counts(oli)
cts_binary = cts>0
sums = rowSums(cts_binary)
det.rate = sums/ncol(cts_binary)
df = as.data.frame(det.rate)

# make plots
print('plotting')
pdf('../plots/Extended_5/distribution_apoe.pdf', width = 3, height = 3)
ggplot(df, aes(x = det.rate)) +
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white", binwidth =0.01 ) +
  geom_density() + geom_vline(xintercept=(c(df['APOE',1])), size=.5, color="black")+theme_classic()
dev.off()

d = as.data.frame(logcounts(oli)['APOE',])
colnames(d) = c('APOE')
clust = kmeans(d$APOE, 2)
clust = (clust$cluster)
d$clust = clust

pdf('../plots/Extended_6/Apoe_density.pdf', width = 3, height = 3)
ggplot(d, aes(x = APOE)) +
  geom_histogram(
                 colour = 1, fill = "white", binwidth =0.1 ) +
  geom_density(aes(fill = clust, alpha = 0.5)) +theme_classic() +facet_wrap( ~ clust, nrow = 3, scales = 'free_y')
dev.off()

# save the info regarding APOE expression/detection
df = df[!df[,1]==0,,drop=F]
df$rank = rank(df[,1])
percentile = df['APOE','rank']/max(df$rank)
out = c(df['APOE',1],percentile*100, table(clust))
names(out) = c('APOE_detection_rate', 'APOE_detection_rate_percentile', 'APOE_detected_pop', 'APOE_not_detected_pop')

write.csv(as.data.frame(out), '../data/supplementary_tables/APOE_expression_post_mortem_oligos.csv')
print('done.')
