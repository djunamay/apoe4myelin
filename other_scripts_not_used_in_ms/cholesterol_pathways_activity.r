library(GSVA)
library(reshape2)

source('../functions/pathway_analyses.r')

pathways = readRDS('../data/other_analyses_outputs/pathways.rds')
low_removed_bp = pathways$pathways$all
cholest_keys =  c('cholest')
cols = readRDS('../data/single_cell_data/Cell_group_colors.rds')

# get genesets that contain the term 'cholesterol'
cholest_gsets = get_gset_by_category(cholest_keys, low_removed_bp)
av_expression = readRDS('../data/single_cell_data/individual_level_averages_per_celltype.rds')

# compute gsva scores for each geneset of the cholesterol type
av = lapply(names(av_expression), function(x) rowMeans(av_expression[[x]]))
names(av) = names(av_expression)
av = do.call('cbind', av)
out = gsva(as.matrix(av), cholest_gsets, mx.diff=TRUE, verbose=TRUE, kcdf=c("Gaussian"), min.sz=5)

df = melt((out))

# plot the distributions of these gsva scores per celltype
df$Var2 = factor(df$Var2, levels = c('Oli', 'Ast', 'Opc','Mic','In', 'Ex'))

#pdf('../figure_outputs/cholesterol_activity.new.pdf', width = 4, height = 5)
boxplot(value ~ Var2, df, horiz = F, las = 3, ylab = 'cholesterol pathways activity', col = unname(cols[levels(df$Var2)]))
#dev.off()

print('done.')
