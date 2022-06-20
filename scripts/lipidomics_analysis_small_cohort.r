# balance the lipidomic cohort according to main covariates

# comparison of APOE3 vs APOE4
library(matrixStats)
# import the human and the ipsc lipidomics data
hdata = read.csv('../external_datasets/small_lipidomics_dataset/Lipidomics_RawData.csv')
hdata2 = read.csv('../external_datasets/small_lipidomics_dataset/Lipidomics_RawData_2.csv')

# combine the separate human lipid data sheets
df = hdata[,c('Molecular.Species','Lipid.Classes',colnames(hdata)[startsWith(colnames(hdata), 'X')])]
df = df[(df$Lipid.Classes)!='',] # remove lipids that don't have an annotation
df2 = hdata2[,c('Molecular.Species','Lipid.Classes',colnames(hdata2)[startsWith(colnames(hdata2), 'X')])]
df2 = df2[(df2$Lipid.Classes)!='',]
out = rbind(df[,colnames(df2)],df2)

# remove zero variance lipids
df = as.data.frame(t(apply(out[,!colnames(out)%in%c('Molecular.Species','Lipid.Classes')], 1, as.numeric)))
colnames(df) = colnames(out)[!colnames(out)%in%c('Molecular.Species','Lipid.Classes')]
keep = !rowVars(as.matrix(df))==0
df = df[keep,]
out = out[keep,]

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

# get differential lipids
apoe = as.data.frame(ifelse((startsWith(as.character(colnames(t)),'4') | startsWith(as.character(colnames(t)),'34')), 'E4', 'E3'))
colnames(apoe) = c('apoe')

# show boxplot of sterol-related changes
ff = t[startsWith(as.character(rownames(t)), 'Sterol'),]
colnames(ff) = paste0(colnames(ff), '_', 1:length(colnames(ff)))

library(reshape2)
library(ggpubr)

#colnames(ff) = c('E4', 'E4','E4','E3','E3','E3')
ff$name = rownames(ff)
fff = melt(ff)
fff$grp = ifelse((startsWith(as.character(fff$variable),'4') | startsWith(as.character(fff$variable),'34')), 'E4', 'E3')
fff$grp = factor(fff$grp, levels = c('E3','E4'))

fff$value = as.numeric(fff$value)
fff$name = as.character(names[fff$name,'name'])
fff = fff[!fff$name=='18:1-d7 Cholesteryl ester',]
fff$name = factor(fff$name, levels = c('20:4 Cholesteryl ester', '18:2 Cholesteryl ester', '22:6 Cholesteryl ester', '18:1 Cholesteryl ester'))
pdf('../plots/human_lipidomics.pdf', width = 4, height = 4)
ggplot(fff, aes(x=grp, y=value, col = grp)) +
      geom_boxplot() + geom_jitter() + facet_grid(. ~ name)+ stat_compare_means(method = "wilcox.test") + theme_classic()
dev.off()
