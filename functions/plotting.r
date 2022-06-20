
####################################################################################################################
#' Filters genesets based on expression: Function only keeps pathways, for which >33% of genes are present in expression vector.
#' For pathways that meet this criterium, genes that are not present in the expression set are removed from that pathway.
#' Finally, pathways with duplicate gene composition (after expression filtering), are removed
#'
#' @param per-celltype list of expressed genes
#' @param list of pathways to filter
#' @return return per-cell-type list of pathways filtered by expression
#' @export

get_go_hmap = function(matrix, names, path_names_go, path_slot){
    temp = as.data.frame((abs(matrix)>1.3) * sign(matrix))
    temp_order = rownames(temp[order((temp$Ex),(temp$In), (temp$Ast), (temp$Mic), (temp$Oli), (temp$Opc), decreasing = T),])
    matrix = matrix[temp_order,]
    index = unname(unlist(lapply(names, function(x) which(rownames(matrix)==x))))
    ha = rowAnnotation(foo = anno_mark(at = index, labels = path_names_go[rownames(matrix),path_slot][index]),gp = gpar(fontsize = 30))
    h = Heatmap(as.matrix(matrix[,c('Ex', 'In', 'Ast', 'Mic', 'Oli', 'Opc')]), right_annotation = ha, col = colorRamp2(c(-2,0,2), c('blue','white','red')),show_row_names = F,border = T, cluster_rows = F,cluster_columns = F,  row_km = 0,rect_gp = gpar(col = 'black', lwd = .5))
    return(h)
}
