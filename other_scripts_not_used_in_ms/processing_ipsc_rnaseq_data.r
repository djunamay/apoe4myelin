########## extended script 1 in APOE4 impairs myelination via cholesterol dysregulation in oligodendrocytes ##########
#############################################################################################################
## iPSC APOE4 vs APOE3 data downloaded from the biomicro servers: 191210Tsa/
# /net/bmc-pub10/data1/bmc/public/Tsai/191210Tsa

source('../functions/bulk_degs.r')

### required libraries
library("edgeR")
library("limma")

# load the counts data
path1 = '../data/iPSC_data/ipsc_bulk_rnaseq_count_files'
file_term = 'geneexp.txt'
counts_no_drug = load_counts(path1, file_term)

# name the counts columns
meta1 = read.table('../data/iPSC_data/ipsc_bulk_rnaseq_count_files/ipsc_metadata.csv', sep = ',', header = F)
rownames(meta1) = meta1$V1

# normalize the counts
norm_counts_no_drug = normalize_counts(counts_no_drug)
write.csv(norm_counts_no_drug, '../data/iPSC_data/ipsc_bulk_no_drug.csv')

# get the DEGs
meta1 = read.table('../data/iPSC_data/ipsc_metadata.csv', sep = ',', header = F)
x = as.data.frame(read.table('../data/iPSC_data/ipsc_bulk_rnaseq_count_files/D19-13224-4367S_geneexp.txt', header = FALSE))
gene.df = as.data.frame(cbind(as.character(x$V1), as.character(x$V3)))
rownames(gene.df) = gene.df$V1
var = as.character(meta1$V2)
var = ifelse(var=='APOE4', 1, 0)
index1 = startsWith(rownames(counts_no_drug), 'ENSG')
counts_no_drug = counts_no_drug[index1,]
index2 = unname(rowSums(counts_no_drug>0)>(0.5*ncol(counts_no_drug)))
counts_no_drug = counts_no_drug[index2,]
out = RunDiffExprAnalysisLimma(counts_no_drug, var, gene.df)
write.csv(out, '../data/iPSC_data/apoe3_v_apoe4_degs.csv')

print('done.')
