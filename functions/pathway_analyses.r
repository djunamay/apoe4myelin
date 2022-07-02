##############################################################################################################################
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
##############################################################################################################################
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

##############################################################################################################################
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

##############################################################################################################################
get_gset_by_category = function(cat, gsets){
  gset = unlist(lapply(names(gsets), function(x) unlist(sum(sapply(cat, grepl, x))>0)))
  gset = (gsets[gset])
  return(gset)
}
