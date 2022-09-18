########## extended script 1 in APOE4 impairs myelination via cholesterol dysregulation in oligodendrocytes ##########
#############################################################################################################
## iPSC APOE4 vs APOE3 data downloaded from the biomicro servers: 191210Tsa/
# /net/bmc-pub10/data1/bmc/public/Tsai/191210Tsa
# TODO: get var info to compute DEGs
# TODO: remove lowly expressed genes

# load_counts function
load_counts = function(path1, file_term){
    files1 = list.files(path1)
    files1 = files1[grepl(file_term, files1)]

    x = read.table(paste0(path1,'/',files1[1]), header = FALSE)
    genenames = x$V1
    genenames = genenames[unlist(lapply(genenames, function(x) startsWith(as.character(x), 'ENSG')))]
    out = list()
    for(i in 1:length(files1)){
      curr = read.table(paste0(path1,'/',files1[i]), header = FALSE)
      rownames(curr) = curr$V1
      curr = curr[as.character(genenames),]
      out[[i]] = curr$V2
    }

    all_counts = do.call('cbind', out)
    rownames(all_counts) = genenames
    colnames(all_counts) = files1

    return(all_counts)
}

# load the counts data
path1 = '../data/iPSC_data/ipsc_bulk_rnaseq_count_files'
file_term = 'geneexp.txt'
counts = load_counts(path1, file_term)
x = strsplit(colnames(counts), '-')
colnames(counts) = unlist(lapply(x, function(i) i[[2]]))

meta1 = read.table('../data/iPSC_data/ipsc_bulk_rnaseq_count_files/ipsc_metadata.csv', sep = ',', header = F)
x = strsplit(meta1$V1, '-')
xx = unlist(lapply(x, function(i) i[[2]]))
x = strsplit(xx, ' :')
rownames(meta1) = unlist(lapply(x, function(i) i[[1]]))
colnames(counts) = meta1[as.character(colnames(counts)),'V2']

# add the gene names
x = as.data.frame(read.table('../data/iPSC_data/ipsc_bulk_rnaseq_count_files/D19-13224-4367S_geneexp.txt', header = FALSE))
gene.df = as.data.frame(cbind(as.character(x$V1), as.character(x$V3)))
rownames(gene.df) = gene.df$V1

index1 = startsWith(rownames(counts), 'ENSG')
counts = as.data.frame(counts[index1,])
counts$gene = as.character(gene.df[rownames(counts), 'V2'])
counts = counts[!duplicated(counts$gene),]
rownames(counts) = counts$gene
counts$gene = NULL
counts = counts[!rowSums(counts)==0,]

write.csv(counts, '../data/iPSC_data/ipsc_opc_rnaseq_counts.txt')