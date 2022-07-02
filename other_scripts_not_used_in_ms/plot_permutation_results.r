
library(ggplot2)

# plot the gset permutation results
scores = readRDS('../data_outputs/pathway_scores.rds')


for(i in unique(names(x))){
  celltype = i
  df = (scores$res$apoe_associated$scores$all[[celltype]])
  observation = nrow(df[df$P.Value<0.05,])

  y = as.data.frame(unlist(x[names(x)==celltype]))
  colnames(y) = 'v1'
  y$rank = rank(y$v1)
  y = y[order(y$v1,decreasing = T),]
  y$percentile = 1-(y$rank/nrow(y))
  fun = approxfun(x = y$v1, y = y$percentile)

  pval = fun(observation)

  plot = ggplot(y, aes(x = v1)) +
    geom_histogram(aes(y = ..density..),
                   colour = 1, fill = "white", binwidth =1 ) +
   geom_vline(xintercept=observation, size=.5, color="black")+theme_classic() + annotate("text", x = 7, y = .3, label = paste0('pval = ', pval)) + ggtitle(celltype)

  pdf(paste0('../plots/gset_permutation_', celltype,'.pdf'), width = 5, height = 5)
  print(plot)
  dev.off()
}


# plot label permutations
x = readRDS('../data_outputs/permutation_analysis_label_perms.rds')
x = do.call('rbind', x)

for(i in unique(colnames(x))){
  celltype = i
  df = (scores$res$GO_BP$scores$all[[celltype]])
  observation = nrow(df[df$P.Value<0.05,])

  y = as.data.frame(x[,i,drop=F])
  colnames(y) = 'v1'
  y$rank = rank(y$v1)
  y = y[order(y$v1,decreasing = T),]
  y$percentile = 1-(y$rank/nrow(y))
  fun = approxfun(x = y$v1, y = y$percentile)

  pval = fun(observation)

  plot = ggplot(y, aes(x = v1)) +
    geom_histogram(aes(y = ..density..),
                   colour = 1, fill = "white", binwidth =1 ) +
   geom_vline(xintercept=observation, size=.5, color="black")+theme_classic() + annotate("text", x = 7, y = .3, label = paste0('pval = ', pval)) + ggtitle(celltype)

  pdf(paste0('../plots/label_permutation_', celltype,'.pdf'), width = 5, height = 5)
  print(plot)
  dev.off()
}
