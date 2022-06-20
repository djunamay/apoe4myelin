
####################################################################################################################
#' Filters genesets based on expression: Function only keeps pathways, for which >33% of genes are present in expression vector.
#' For pathways that meet this criterium, genes that are not present in the expression set are removed from that pathway.
#' Finally, pathways with duplicate gene composition (after expression filtering), are removed
#'
#' @param per-celltype list of expressed genes
#' @param list of pathways to filter
#' @return return per-cell-type list of pathways filtered by expression
#' @export

filter_lowly_exp_genes = function(expressed, all_paths){
    gsets_per_celltype = list()

    for(i in names(expressed)){
        print(i)
        index = unlist(lapply(names(all_paths), function(x) (sum(all_paths[[x]]%in%expressed[[i]])/length(all_paths[[x]]))>0.33))
        p = all_paths[(index)]
        x = lapply(names(p), function(y) intersect(expressed[[i]], p[[y]]))
        names(x) = names(p)
        x = x[!duplicated(x)]
        gsets_per_celltype[[i]] = x

    }
    return(gsets_per_celltype)
}

####################################################################################################################
#' get list of marker genes
#'
#' @param species
#' @param category
#' @param subcategory
#' @return return GO terms
#' @export
get_fits = function(gsva_out, meta){
    fits = list()
    for(i in names(gsva_out)){
        predict = meta[as.character(rownames(gsva_out[[i]])),]
        mod = model.matrix(~APOE4 + amyloid + nft + age_death + msex + pmi, data=predict)
        fits[[i]] = fit.gsva(mod, i, gsva_out, 'APOE4e4')

    }
    return(fits)
}

####################################################################################################################
#' get list of marker genes
#'
#' @param species
#' @param category
#' @param subcategory
#' @return return GO terms
#' @export
get_scores = function(fits){
    outs = list()
    all = list()
    for(i in names(fits)){
        df = fits[[i]]
        df$score = sign(df$logFC) * -log10(df$P.Value)
        df$celltype = i
        df = df[order(abs(df$logFC),decreasing = T),]
        all[[i]] = df[,c('celltype', 'logFC','P.Value', 'names')]
        df = df[df$P.Value<0.05,]
        if(nrow(df)<10){
            outs[[i]] =  df
        }else{
            outs[[i]] =  df[1:10,c('celltype', 'score', 'names')]
        }
    }
    return(list('all' = all, 'top' = outs))
}

####################################################################################################################
#' get list of marker genes
#'
#' @param species
#' @param category
#' @param subcategory
#' @return return GO terms
#' @export
get_matrix = function(scores, top_paths){
    df = do.call('rbind',scores)
    df$score = sign(df$logFC) * -log10(df$P.Value)
    df = as.data.frame(pivot_wider(df[,c('celltype', 'names', 'score')], values_from = 'score', names_from = 'celltype'))
    df[is.na(df)] = 0
    rownames(df) = df$names
    df$names = NULL
    return(df[top_paths,])
}

# TODO clean up below functions
####################################################################################################################
####################################################################################################################
####################################################################################################################
#' get list of marker genes
#'
#' @param species
#' @param category
#' @param subcategory
#' @return return GO terms
#' @export
get_go_terms = function(species, cat, subcat){
    h_gene_sets = msigdbr::msigdbr(species = species, category = cat, subcategory = subcat)
    gs = list()
    for(i in unique(h_gene_sets$gs_name)){
        temp = h_gene_sets[h_gene_sets$gs_name==i,c('human_gene_symbol','gs_name')]
        #name = c(name,temp$gs_name[1])
        d = unique(temp$human_gene_symbol)
        gs[[i]] = d
    }
    return(gs)
}

####################################################################################################################
#' get list of marker genes
#'
#' @param colors
#' @param celltype
#' @param pval_cut
#' @param pval_slot
#' @return volcano plot
#' @export
get_hmap_pathways = function(df, colors, celltype_vector, pval_cut, pval_slot, pathway_slot, NES_slot, names){
    cols = c('grey',colors[celltype_vector])
    names(cols) = c('other',celltype_vector)
    c = ifelse(df[[pval_slot]]<pval_cut, celltype_vector, 'other')

    if(names==TRUE){
        plot = ggplot2::ggplot(data = df, ggplot2::aes(x = .data[[NES_slot]], y = -log10(.data[[pval_slot]]))) +
        ggplot2::geom_point(ggplot2::aes( col = c),size=5) + ggplot2::theme_classic() +
        ggrepel::geom_text_repel(label = ifelse(df[[pval_slot]]<pval_cut, as.data.frame(df)[,pathway_slot], ''), size= 2,force=5, min.segment.length=1,force_pull= 0, max.overlaps = 30) +
        ggplot2::scale_color_manual(values = cols) + ggplot2::theme(legend.position="none")  #+ geom_hline(yintercept=3, linetype="dashed", color = "black")
        return(plot)
    }else{
        plot = ggplot2::ggplot(data = df, ggplot2::aes(x = .data[[NES_slot]], y = -log10(.data[[pval_slot]]))) +
        ggplot2::geom_point(ggplot2::aes( col = c),size=5) + ggplot2::theme_classic() +
        ggplot2::scale_color_manual(values = cols) + ggplot2::theme(legend.position="none")   #+ geom_hline(yintercept=3, linetype="dashed", color = "black")
        return(plot)
    }

}

####################################################################################################################
check_rownames = function(h.gsva){
    prev = rownames(h.gsva[[1]])
    out = c()
    for (i in names(h.gsva)){
      curr = rownames(h.gsva[[i]])
      #print(unique(prev == curr))
      out = c(out, sum(curr == prev))
     }
    if(sum(out)==length(h.gsva)*length(prev)){
        print('the sample orders are identical')
    }else{
        print('the sample orders are different')
    }

}
####################################################################################################################
 fit.gsva = function(mod1, i, gsva.per.celltype, coef){
  # fit linear model to gsva pathway activity scores
    fit <- lmFit(t(gsva.per.celltype[[i]]), design=mod1)
    fit <- eBayes(fit)
    allgenesets <- topTable(fit, coef=coef, number=Inf, confint = T) %>% .[order(.$P.Value, decreasing = F),]
    allgenesets$celltype = i
    allgenesets$names = rownames(allgenesets)
    return(allgenesets)

  }
####################################################################################################################
get_gset_names_by_category = function(cat, gsets){


  gset = unlist(lapply(gsets, function(x) unlist(sum(sapply(cat, grepl, x))>0)))
  gset = (gsets[gset])

  return(gset)

}
####################################################################################################################
#' compute hypergeometric p-values on geneset overlap with signature
#'
#' @param signature: vector of genes for which to compute enrichment
#' @param genesets: list of gene sets
#' @param background: length of vector (genes) from which signature was sampled
#' @return hyp enrichment
#' @export
hyp = function(signature, genesets, background){

      signature.found <- signature[signature %in% unique(unlist(genesets))]
      n.hits <- sapply(genesets, function(x, y) length(intersect(x, y)), signature.found)
      n.drawn <- length(signature)
      n.genesets <- sapply(genesets, length)
      n.left <- background-n.genesets

      pvals <- suppressWarnings(stats::phyper(q=n.hits-1,
                                            m=n.genesets,
                                            n=n.left,
                                            k=n.drawn,
                                            lower.tail=FALSE))

      data <- data.frame(label=names(genesets),
                       pval=signif(pvals, 2),
                       fdr=signif(stats::p.adjust(pvals, method="fdr"), 2),
                       signature=length(signature),
                       geneset=n.genesets,
                       overlap=n.hits,
                       background=background,
                       hits=sapply(genesets, function(x, y) paste(intersect(x, y), collapse=','), signature.found),
                       stringsAsFactors=FALSE)
    return(data)
}

get_sem_sim_matrix = function(go_terms){
    # input: vector of go names
    all_combs = t(combn(go_terms, 2))
    all_combs = as.data.frame(all_combs)
    all_combs2 = cbind(as.character(all_combs$V2), as.character(all_combs$V1))
    colnames(all_combs2) = c('V1', 'V2')
    all_combs = rbind(all_combs, all_combs2)

    all_combs$sim = unlist(lapply(1:nrow(all_combs), function(x) goSim(as.character(all_combs[x,1]),as.character(all_combs[x,2]), semData=hsGO)))
    #all_combs = as.data.frame(all_combs)
    #all_combs$V1 = as.character(all_combs$V1)
    #all_combs$V2 = as.character(all_combs$V2)


    #df = edgelist_to_adjmat(all_combs[,1:2], w = as.numeric(as.character(all_combs[,3])), t0 = NULL,t1 = NULL,t = NULL, simplify = TRUE, undirected = TRUE, self = TRUE, multiple = TRUE, keep.isolates = TRUE, recode.ids = TRUE)
    #df = as.matrix(df)[go_terms,go_terms]
    return(all_combs)
}





# show visualization of the top pathways


get_gene_matrix = function(path_names, genesets){
    o = list()
    d = genesets[path_names]
    for(i in names(d)){
        x = as.data.frame(d[i])
        colnames(x) = c('genes')
        x$path = i
        x$present = 1
        o[[i]] = x
    }
    o = do.call('rbind',o)
    d = as.data.frame(pivot_wider(o, values_from = 'present', names_from = 'path'))
    rownames(d) = d$genes
    d$genes = NULL
    d[is.na(d)] = 0
    return(d)
}
