##############################################################################################################################
get_pathway_fits = function(order, av_expression, pathways, top_20, summary){
    # run GSVA on GO pathways
    all_data = list()
    print('running GSVA...')
    out_bp_terms = lapply(order, function(x) t(gsva(as.matrix(av_expression[[x]]), pathways[[x]], mx.diff=TRUE, verbose=F, kcdf=c("Gaussian"), min.sz=5, max.sz = 150, parallel.sz=16)))
    names(out_bp_terms) = order
    all_data[['gsva_out']] = out_bp_terms

    # get linear model fits
    print('getting linear model fits...')
    fits = get_fits(out_bp_terms, summary)
    all_data[['fits_all']] = fits

    # get matrix of scores for heatmap
    print('get matrix of scores')
    scores = get_scores(fits)
    all_data[['scores_all']] = scores

    print('filter by score 1.3')
    names = unique(unname(unlist(lapply(names(scores$all), function(x) rownames(scores$all[[x]])))))
    mat = get_matrix(scores$all, names)
    mat = mat[unname(rowSums(abs(mat)>1.3)>0),]

    if(top_20==TRUE){
        index = unique(unname(unlist(lapply(colnames(mat), function(x) order(abs(mat[[x]]),decreasing = T)[1:20]))))
        mat = mat[index,]
    }
    all_data[['scores_filtered']] = mat
    return(all_data)
}
####################################################################################################################
get_fits = function(gsva_out, meta){
    fits = list()
    for(i in names(gsva_out)){
        predict = meta[as.character(rownames(gsva_out[[i]])),]
        mod = model.matrix(~APOE4 + amyloid + nft + age_death + msex + pmi, data=predict)
        fits[[i]] = fit.gsva(mod, i, gsva_out, 'APOE4')
    }
    return(fits)
}
####################################################################################################################
 fit.gsva = function(mod1, i, gsva.per.celltype, coef){
    fit <- lmFit(t(gsva.per.celltype[[i]]), design=mod1)
    fit <- eBayes(fit)
    allgenesets <- topTable(fit, coef=coef, number=Inf, confint = T) %>% .[order(.$P.Value, decreasing = F),]
    allgenesets$celltype = i
    allgenesets$names = rownames(allgenesets)
    return(allgenesets)

  }
####################################################################################################################
get_scores = function(fits){
    outs = list()
    all = list()
    for(i in names(fits)){
        df = fits[[i]]
        df$celltype = i
        df = df[order(abs(df$logFC),decreasing = T),]
        all[[i]] = df[,c('celltype', 'logFC','P.Value', 'names')]
    }
    return(list('all' = all))
}
####################################################################################################################
get_matrix = function(scores, top_paths){
    df = do.call('rbind',scores)
    df$score = sign(df$logFC) * -log10(df$P.Value)
    df = as.data.frame(pivot_wider(df[,c('celltype', 'names', 'score')], values_from = 'score', names_from = 'celltype'))
    df[is.na(df)] = 0
    rownames(df) = df$names
    df$names = NULL
    return(df[top_paths,])
}
####################################################################################################################
get_stratified_fits = function(stratified_summary, gsva_out_full, model, coefficient){
    out = lapply(names(gsva_out_full), function(x) gsva_out_full[[x]][rownames(gsva_out_full[[x]])%in%as.character(rownames(stratified_summary)),])
    names(out) = names(gsva_out_full)
    fit <- lmFit(t(out[['Oli']]), design=model)
    fit <- eBayes(fit)
    allgenesets <- topTable(fit, coef=coefficient, number=Inf, confint = T) %>% .[order(.$P.Value, decreasing = F),]
    return(allgenesets)
}
