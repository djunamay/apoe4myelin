##############################################################################################################################
filter_lowly_exp_genes = function(expressed, all_paths){
    gsets_per_celltype = list()

    for(i in names(expressed)){
        index = unlist(lapply(names(all_paths), function(x) (sum(all_paths[[x]]%in%expressed[[i]])/length(all_paths[[x]]))>0.33))
        p = all_paths[(index)]
        x = lapply(names(p), function(y) intersect(expressed[[i]], p[[y]]))
        names(x) = names(p)
        x = x[!duplicated(x)]
        gsets_per_celltype[[i]] = x
    }
    return(gsets_per_celltype)
}
##############################################################################################################################
read.geneset = function(path_to_gset){
  bp = GSA.read.gmt(path_to_gset)
  out = bp$genesets
  out = lapply(1:length(out), function(x) out[[x]][out[[x]]!=''])
  names(out) = bp$geneset.names
  return(out)
}
####################################################################################################################
get_gset_names_by_category = function(cat, gsets){
  gset = unlist(lapply(gsets, function(x) unlist(sum(sapply(cat, grepl, x))>0)))
  gset = (gsets[gset])
  return(gset)
}
####################################################################################################################
