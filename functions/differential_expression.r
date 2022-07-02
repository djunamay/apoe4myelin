##############################################################################################################################
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
##############################################################################################################################
RunDiffExprAnalysisLimma <- function(counts.df, var, genes.df) {
    # Computes differential expression using a combination of sva and voom-limma
    #
    # Args:
    #   counts.df: data.frame with read counts (#genes x #samples)
    #   var: variable of interest
    #   genes.df: data.frame mapping gene IDs to gene names
    #   n.sv: number of surrogates for SVA. Automatically determined if NULL.
    #
    # Returns:
    #   voom-limma output augmented with Storey q-values
    #   SVA surrogate variables

    design <- cbind(1, var)

    # apply edgeR normalization (TMM) to counts
    dge <- DGEList(counts=counts.df) # makes a DGEList object from the table of counts
    dge <- calcNormFactors(dge) # adjusts for RNA composition (adjusting for differences in per cell total RNA output)

    # re-calculate voom transformation w/ SVA covariates; run limma
    mod1 <- model.matrix(~var)  # intercept, var, SVs
    v <- voom(dge, design=mod1)
    fit <- lmFit(v, design=mod1)
    fit <- eBayes(fit)
    res <- topTable(fit, coef=ncol(design), n=Inf, sort.by="none")

    res[, 'gene_name'] <- genes.df[row.names(res), 'V2']
    #res <- res[, colnames(res)[c(8,1:7)]]

    return(list("res"=res))
}
##############################################################################################################################
WriteResults <- function(data, filename, header) {
    datafile <- file(filename, open = 'wt')
    on.exit(close(datafile))
    writeLines(paste0(header, collapse="\t"), con=datafile, sep='\n')
    write.table(data, datafile, sep='\t', col.names=FALSE, quote=FALSE)
}
##############################################################################################################################
normalize_counts = function(counts.df){
  # apply edgeR normalization (TMM) to counts
  dge <- DGEList(counts=counts.df) # makes a DGEList object from the table of counts
  dge <- calcNormFactors(dge) # adjusts for RNA composition (adjusting for differences in per cell total RNA output)
  normalized = cpm(dge)
  return(normalized)
}
##############################################################################################################################
Wilcox.differential <- function(QInexp, QIngroup) {
 o <- presto::wilcoxauc(QInexp, QIngroup)
 o <- o[o$group==1,]
 o <- o[order(o$auc, decreasing = T),]
 rownames(o) <- o$feature
 o$class.power <- (abs(o$auc-0.5)) * 2
 return(o)
}
##############################################################################################################################
RunDiffExprAnalysisLimma_pseudobulk <- function(logcounts, var, predict) {
    mod1 = model.matrix(~APOE4 + amyloid + nft + age_death + msex + pmi, data=predict)
    fit <- lmFit(logcounts, design=mod1)
    fit <- eBayes(fit)
    res <- topTable(fit, coef=var, n=Inf, sort.by="none")
    return(list("res"=res))
}
##############################################################################################################################
lipid_species_wilcox_test = function(data_subset,grp){
    stats = list()
    for(i in colnames(data_subset)[2:5]){
      y = data_subset[data_subset$apoe_genotype%in%c('33','23','22'),i]
      x = data_subset[!data_subset$apoe_genotype%in%c('33','23','22'),i]
      temp = wilcox.test(x,y,conf.int=T, conf.level = 0.90)
      stats[[i]] = c(temp$conf.int, temp$estimate, temp$p.value)
    }
    stats = do.call('rbind', stats)
    stats = as.data.frame(stats)
    colnames(stats) = c('lower_90_CI', 'upper_90_CI', 'sample_estimates.median_difference_APOE4carrier_vs_noncarrier.', 'p.value')
    stats$grp = grp
    return(stats)
}
