##############################################################################################################################
get_go_hmap = function(matrix, names, path_names_go, path_slot){
    temp = as.data.frame((abs(matrix)>1.3) * sign(matrix))
    temp_order = rownames(temp[order((temp$Ex),(temp$In), (temp$Ast), (temp$Mic), (temp$Oli), (temp$Opc), decreasing = T),])
    matrix = matrix[temp_order,]
    index = unname(unlist(lapply(names, function(x) which(rownames(matrix)==x))))
    ha = rowAnnotation(foo = anno_mark(at = index, labels = path_names_go[rownames(matrix),path_slot][index]),gp = gpar(fontsize = 30))
    h = Heatmap(as.matrix(matrix[,c('Ex', 'In', 'Ast', 'Mic', 'Oli', 'Opc')]), right_annotation = ha, col = colorRamp2(c(-2,0,2), c('blue','white','red')),show_row_names = F,border = T, cluster_rows = F,cluster_columns = F,  row_km = 0,rect_gp = gpar(col = 'black', lwd = .5))
    return(h)
}
##############################################################################################################################
boxplot_w_stats = function(df, x, y, group_color = x, group_fill = x, alpha=.5, palette, xlab='', ylab='', width=.5, stats_method = 'wilcox', comparisons){
    plot <- ggpubr::ggboxplot(df, x = x, y = y,
          color = group_color, fill = group_fill, alpha = alpha, palette = palette,
          xlab = xlab, ylab = ylab, width = width)+ ggpubr::stat_compare_means(method = stats_method, comparisons = comparisons) + ggplot2::geom_point(alpha = .3) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))
    return(plot)

}
##############################################################################################################################
get_barplot = function(all_data, x, y){

    df = as.matrix(t(table(all_data[,x], all_data[,y])))
    names = unique(all_data[,y])
    coul <- RColorBrewer::brewer.pal(length(names), "Set2")
    names(coul) = names
    res = stats::fisher.test(df)

    bp = graphics::barplot(t(apply(df, 1, function(x){x/colSums(df)})), col = coul, legend = TRUE, las = 2)
    graphics::text(bp, paste0('p-value = ',round(res$p.value, digits = 2)), ylab = y)
    return(bp)

}
