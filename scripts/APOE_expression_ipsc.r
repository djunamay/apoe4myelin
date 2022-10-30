########## APOE expression related to extended data figure 6 #############
###################################################################
print('|| plotting APOE expression in oligodendrocytes (iPSC)... ||')

# required libraries
library(ggplot2)

counts = read.csv('../data/iPSC_data/FPKM_table_OPC.txt', sep = '\t')
rownames(counts) = counts$gene
counts$gene = NULL

# get percentile in which APOE is expressed
means = rowMeans(counts)
df = as.data.frame(means)
df = df[!df$means==0,,drop = F]

d = df[order(df$means,decreasing = F),,drop = F]
percentile = (which(rownames(d) =='APOE')/nrow(d))*100

pdf('../plots/Extended_5/apoe_expression_ipsc.pdf', width = 5, height = 5)
ggplot(df, aes(x = (means))) +
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white", binwidth =0.5 ) +
  geom_density() + geom_vline(xintercept=(c((df['APOE',1]))), size=.5, color="black")+theme_classic() + xlim(c(0,50)) + ylim(c(0,.27))+ annotate("text", x = 10, y=0.25, label = percentile)
dev.off()

print('done.')
