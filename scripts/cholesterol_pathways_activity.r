# get genesets that contain the term 'cholesterol'
cholest_gsets = get_gset_by_category(cholest_keys, all)

# compute gsva scores for each geneset of the cholesterol type
av = lapply(names(av_expression), function(x) rowMeans(av_expression[[x]]))
names(av) = names(av_expression)
av = do.call('cbind', av)
out = gsva(as.matrix(av), cholest_gsets, mx.diff=TRUE, verbose=TRUE, kcdf=c("Gaussian"), min.sz=5)
df = melt(out)

# plot the distributions of these gsva scores per celltype
options(repr.plot.height=5, repr.plot.width=4)
df$X2 = factor(df$X2, levels = c('Oli', 'Ast', 'Mic','Opc','In', 'Ex'))

#pdf('../figure_outputs/cholesterol_activity.new.pdf', width = 4, height = 5)
boxplot(value ~ X2, df, horiz = F, las = 3, ylab = 'cholesterol pathways activity', col = unname(cols[levels(df$X2)]))
#dev.off()
