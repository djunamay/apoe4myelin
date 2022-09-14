########## analysis of lipidomic data from human corpus callosum ##########
###########################################################################

print('|| analysis of post-mortem CC lipidomic data... ||')

# required packages
library(reshape2)
library(ggpubr)
library(matrixStats)

# import the human post-mortem lipidomics data
print('loading and processing the data')
hdata = read.csv('../data/lipidomics_dataset/cc_lipidomics/Lipidomics_RawData.csv')
hdata2 = read.csv('../data/lipidomics_dataset/cc_lipidomics/Lipidomics_RawData_2.csv')

# combine the separate human lipid data sheets
df = hdata[,c('Molecular.Species','Lipid.Classes',colnames(hdata)[startsWith(colnames(hdata), 'X')])]
df = df[(df$Lipid.Classes)!='',] # remove lipids that don't have an annotation
df2 = hdata2[,c('Molecular.Species','Lipid.Classes',colnames(hdata2)[startsWith(colnames(hdata2), 'X')])]
df2 = df2[(df2$Lipid.Classes)!='',]
out = rbind(df[,colnames(df2)],df2)

df = as.data.frame(t(apply(out[,!colnames(out)%in%c('Molecular.Species','Lipid.Classes')], 1, as.numeric)))
colnames(df) = colnames(out)[!colnames(out)%in%c('Molecular.Species','Lipid.Classes')]

# extract the experimental group
index = colnames(df)
index = strsplit(index, '[.]')
grp = unlist(lapply(1:length(index), function(x) tail(index[[x]], n = 1)))

# set the column and row names
colnames(df) = grp
rownames(df) = paste0(out$Lipid.Classes,'_', 1:length(out$Lipid.Classes))
names = cbind(as.data.frame(out[,'Molecular.Species']),as.data.frame(rownames(df)))
colnames(names) = c('name', 'id')
rownames(names) = names$id

# subset by sex
t = df[,endsWith(colnames(df), 'F')]

# remove zero variance lipids
keep = !rowVars(as.matrix(t))==0
t = t[keep,]
names = names[keep,]

# get differential lipids
apoe = as.data.frame(ifelse((startsWith(as.character(colnames(t)),'4') | startsWith(as.character(colnames(t)),'34')), 'E4', 'E3'))
colnames(apoe) = c('apoe')

# prep data for boxplot
print('plotting the data')
ff = t[startsWith(as.character(rownames(t)), 'Sterol'),]
colnames(ff) = paste0(colnames(ff), '_', 1:length(colnames(ff)))
ff$name = rownames(ff)
fff = melt(ff)
fff$grp = ifelse((startsWith(as.character(fff$variable),'4') | startsWith(as.character(fff$variable),'34')), 'E4', 'E3')
fff$grp = factor(fff$grp, levels = c('E3','E4'))
fff$value = as.numeric(fff$value)
fff$name = as.character(names[fff$name,'name'])
fff = fff[!fff$name=='18:1-d7 Cholesteryl ester',]
fff$name = factor(fff$name, levels = c('20:4 Cholesteryl ester', '18:2 Cholesteryl ester', '22:6 Cholesteryl ester', '18:1 Cholesteryl ester'))

# draw boxplot
pdf('../plots/Figure_2/human_lipidomics_cc.pdf', width = 4, height = 4)
ggplot(fff, aes(x=grp, y=value, col = grp)) +
      geom_boxplot() + geom_jitter() + facet_grid(. ~ name)+ stat_compare_means(method = "wilcox.test") + theme_classic()
dev.off()

# save wilcox p-values and effect sizes for each lipid species
stats = list()
for(i in 1:nrow(t)){
  x = unlist(t[i,1:3])
  y = unlist(t[i,4:6])
  temp = wilcox.test(x,y,conf.int=T, conf.level = 0.90)
  stats[[i]] = c(temp$conf.int, temp$estimate, temp$p.value)
}

stats = do.call('rbind', stats)
stats = as.data.frame(stats)
colnames(stats) = c('lower_90_CI', 'upper_90_CI', 'sample_estimates(median_difference_44_vs_33)', 'p.value')
rownames(stats) = rownames(t)
stats$names = names[rownames(stats), 'name']
all_data = cbind(t, stats)
write.csv(all_data, '../data/supplementary_tables/cc_lipidomics_data.csv')
print('done.')
