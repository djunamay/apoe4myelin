#################################################################################
# auxiliary functions
#################################################################################
QC.filter.sce <- function(sce, nGenesL=500, nGenesH=10000, nCellsObs=10) {
  print("initial dim")
  print(dim(sce))
  nGenesDetected <- Matrix::colSums(assays(sce)[["counts"]]>0)
  nCellsDetected <- Matrix::rowSums(assays(sce)[["counts"]]>0)
  sce@colData$nGenesDetected <- nGenesDetected
  rowData(sce)$nCellsDetected <- nCellsDetected
  which(sce@colData$nGenesDetected>=nGenesL & sce@colData$nGenesDetected<=nGenesH) -> keep
  ######################################################
  print("cells depth")
  print(c(ncol(sce),length(keep)))
  ######################################################
  sce <- sce[,keep]
  gc()
  rowData(sce)$nCellsDetected <- Matrix::rowSums(assays(sce)[["counts"]]>0)
  
  Genes.complete <- readRDS("../data/single_cell_data/ensembl.GRCh38p12.genes.complete.annot.rds")
  Genes <- Genes.complete[Genes.complete$Gene.type=="protein_coding"]
  Genes <- unique(Genes$Gene.name)
  
  sce <- sce[as.character(rowData(sce)[,2])%in%Genes,]
  gc()
  ######################################################
  print("genes non detected")
  print(sum(rowData(sce)$nCellsDetected<=nCellsObs))
  ######################################################
  sce <- sce[rowData(sce)$nCellsDetected>nCellsObs, ]
  gc()
  sce <- sce[!duplicated(as.character(rowData(sce)[,2])),]
  gc()
  rownames(sce) <- as.character(rowData(sce)[,2])
  rownames(assays(sce)[["counts"]]) <- as.character(rowData(sce)[,2])
  
  Mito <- rownames(sce)[grep("MT-", rownames(sce))][-1]
  TotalCounts <- colSums(assays(sce)[["counts"]])
  MitoCounts <- colSums(assays(sce)[["counts"]][Mito, ])
  NonMitoCounts <- colSums(assays(sce)[["counts"]][!rownames(sce)%in%Mito, ])
  
  sce@colData$TotalCounts <- TotalCounts
  sce@colData$MitoCounts <- MitoCounts
  sce@colData$NonMitoCounts <- NonMitoCounts
  sce@colData$MitoNonMitoRation <- sce@colData$MitoCounts/sce@colData$NonMitoCounts
  sce@colData$MitoNonMitoRationPass <- !Binarize(sce@colData$MitoNonMitoRation)>0
  ######################################################
  print("cells mito")
  print(sum(!sce@colData$MitoNonMitoRationPass))
  ######################################################
  sce <- sce[!rownames(sce)%in%Mito, sce@colData$MitoNonMitoRationPass]
  gc()
  ######################################################
  print("final dim")
  print(dim(sce))
  ######################################################
  return(sce)
}
######################################
Binarize <- function(InScores) {
  X <- kmeans(InScores, 2)
  out <- X$cluster*0
  out[X$cluster==which.max(as.numeric(X$centers))] <- 1
  return(out)
}
######################################
ACTIONet.clustering.analysis <- function(xsce, Renorm=FALSE, doCluster = TRUE, ClustLabel="cluster", ...) {
  if(Renorm) {
    assays(xsce)[[2]] <- NULL
    reducedDim(xsce, type = "S_r") <- NULL
    #reducedDim(xsce)[[1]] <- NULL
    xsce <- reduce.sce(xsce)
  }
  out <- run.ACTIONet(xsce, k_max =2, ...)
  if(doCluster) out <- cluster.ACTIONet(out, annotation.name = ClustLabel)
  
  xsce@metadata$network <- out$ACTIONet
  V(xsce@metadata$network)$name <- colnames(xsce)
  xsce@metadata$vis.out <- out$vis.out
  xsce@metadata[[ClustLabel]] <- paste0("C", as.character(out$annotations[[ClustLabel]]$Labels))
  
  xsce@metadata[[paste0(ClustLabel,".profiles")]] <- Get.column.group.average(x = assays(xsce)[["logcounts"]], group = xsce@metadata[[ClustLabel]])
  xsce@metadata[[paste0(ClustLabel,".signatures")]] <- Aggregate.pairwise.FC.colPairwise(xsce@metadata[[paste0(ClustLabel,".profiles")]])
  xsce@metadata[[paste0(ClustLabel,".top.genes")]] <- Extract.tops.by.column(xsce@metadata[[paste0(ClustLabel,".signatures")]])
  xsce@metadata[[paste0(ClustLabel,".cor")]] <- cor(xsce@metadata[[paste0(ClustLabel,".signatures")]])
  
  DF <- DataFrame(num=1:ncol(xsce), tag=colnames(xsce),
                  x=xsce@metadata$vis.out$coordinates[,1], y=xsce@metadata$vis.out$coordinates[,2],
                  X=xsce@metadata$vis.out$coordinates_3D[,1], Y=xsce@metadata$vis.out$coordinates_3D[,2], Z=xsce@metadata$vis.out$coordinates_3D[,3],
                  color=rgb(xsce@metadata$vis.out$colors),
                  connectivity=V(xsce@metadata$network)$connectivity/max(V(xsce@metadata$network)$connectivity),
                  cluster=xsce@metadata$cluster
  )
  xsce@colData <- DataFrame(xsce@colData, DF)
  rm(out)
  gc()
  return(xsce)
}
######################################
Add.best.annotation <- function(xsce, geneSet, Annotlabel, ClustLabel="cluster") {
  xsce@metadata[[paste0(Annotlabel, ".enrichment")]] <- Permutation.enrichment.analysis(t(xsce@metadata[[paste0(ClustLabel,".signatures")]]), marker.genes = geneSet)
  xsce@colData[,Annotlabel] <- xsce@metadata[[paste0(Annotlabel, ".enrichment")]]$annotation[xsce@colData[,ClustLabel], "labels"]
  return(xsce)
}
######################################
Permutation.enrichment.analysis <- function(x, marker.genes, rand.sample.no = 1000) {
  #require(ACTIONet)
  require(igraph)
  require(Matrix)
  require(stringr)
  
  if(is.matrix(marker.genes) | is.sparseMatrix(marker.genes)) {
    marker.genes = apply(marker.genes, 2, function(x) rownames(marker.genes)[x > 0])
  }
  
  archetype.panel <- x
  
  GS.names = names(marker.genes)
  if (is.null(GS.names)) {
    GS.names = sapply(1:length(GS.names), function(i) sprintf("Celltype %s", i))
  }
  
  markers.table = do.call(rbind, lapply(names(marker.genes), function(celltype) {
    genes = marker.genes[[celltype]]
    if (length(genes) == 0)
      return(data.frame())
    
    
    signed.count = sum(sapply(genes, function(gene) grepl("\\+$|-$", gene)))
    is.signed = signed.count > 0
    
    if (!is.signed) {
      df = data.frame(Gene = (genes), Direction = +1, Celltype = celltype)
    } else {
      
      pos.genes = (as.character(sapply(genes[grepl("+", genes, fixed = TRUE)], function(gene) stringr::str_replace(gene,
                                                                                                                   stringr::fixed("+"), ""))))
      neg.genes = (as.character(sapply(genes[grepl("-", genes, fixed = TRUE)], function(gene) stringr::str_replace(gene,
                                                                                                                   stringr::fixed("-"), ""))))
      
      df = data.frame(Gene = c(pos.genes, neg.genes), Direction = c(rep(+1, length(pos.genes)), rep(-1, length(neg.genes))),
                      Celltype = celltype)
    }
  }))
  markers.table = markers.table[markers.table$Gene %in% colnames(archetype.panel), ]
  
  if (dim(markers.table)[1] == 0) {
    print("No markers are left")
    return()
  }
  
  IDX = split(1:dim(markers.table)[1], markers.table$Celltype)
  
  print("Computing significance scores")
  set.seed(0)
  Z = sapply(IDX, function(idx) {
    markers = (as.character(markers.table$Gene[idx]))
    directions = markers.table$Direction[idx]
    mask = markers %in% colnames(archetype.panel)
    
    A = as.matrix(archetype.panel[, markers[mask]])
    sgn = as.numeric(directions[mask])
    stat = A %*% sgn
    
    rand.stats = sapply(1:rand.sample.no, function(i) {
      rand.samples = sample.int(dim(archetype.panel)[2], sum(mask))
      rand.A = as.matrix(archetype.panel[, rand.samples])
      rand.stat = rand.A %*% sgn
    })
    
    cell.zscores = as.numeric((stat - apply(rand.stats, 1, mean))/apply(rand.stats, 1, sd))
    
    return(cell.zscores)
  })
  
  Z[is.na(Z)] = 0
  Labels = colnames(Z)[apply(Z, 1, which.max)]
  Labels.conf = apply(Z, 1, max)
  
  names(Labels) = rownames(archetype.panel)
  names(Labels.conf) = rownames(archetype.panel)
  rownames(Z) = rownames(archetype.panel)
  
  out.list = list(annotation=data.frame(labels = Labels, labels.confidence = Labels.conf, stringsAsFactors = FALSE), enrichment = Z)
  
  return(out.list)
}
######################################
Extract.enrichments.from.sce <- function(xsce) {
  out <- xsce@metadata[grep(".enrichment", names(xsce@metadata))]
  return(out)
}
######################################
Extract.summarized.pattern.output.from.sce <- function(xsce, metatoExtract=c("cluster.profiles", "cluster.signatures", "cluster.top.genes", "cluster.cor"), getEnrichments=FALSE) {
  out <- xsce@metadata[metatoExtract]
  if(getEnrichments) out$enrichments <- Extract.enrichments.from.sce(xsce)
  out$colMetadata <- xsce@colData
  return(out)
}
######################################
Sub.sce.cols <- function(InSCE, xcIndx) {
  keep <- xcIndx
  InSCE@metadata$network <- induced.subgraph(InSCE@metadata$network, V(InSCE@metadata$network)$name[keep])
  InSCE@metadata$cluster <- InSCE@metadata$cluster[keep]
  InSCE@metadata$vis.out$coordinates <- InSCE@metadata$vis.out$coordinates[keep,]
  InSCE@metadata$vis.out$coordinates_3D <- InSCE@metadata$vis.out$coordinates_3D[keep,]
  InSCE@metadata$vis.out$colors <- InSCE@metadata$vis.out$colors[keep,]
  return(InSCE[,keep])
}
######################################
# plotting functions
Plot.coords <- function(xdf, x="x", y="y", Col="tomato", Trans=1, cex = 0.15, add.text=F, add.legend=F, legCol="black", legbg="gold", legcex=0.8, mapLegColor=TRUE, ...) {
  plot(xdf[,x], xdf[,y], pch = 21, cex = cex, bg = ggplot2::alpha(Col, Trans), col = colorspace::darken(ggplot2::alpha(Col, Trans), 0.25), axes = F, xlab = "", ylab = "", ...)
  
  if ( add.text == TRUE) {
    par(xpd = T, mar = par()$mar * c(1.1,1.1,1.1,1.1))
    
    colMapping <- Col[unique(names(Col))]
    Annot <- names(colMapping)
    labels <- names(Col)
    
    require(wordcloud)
    centroids = t(sapply(Annot, function(l) {
      idx = which(labels == l)
      if(length(idx) == 1) {
        return(as.numeric(xdf[idx, c(x,y)]))
      }
      
      sub.coors = xdf[idx, c(x,y)]
      anchor.coor = as.numeric(apply(sub.coors, 2, function(x) mean(x, trim = 0.80)))
      
      return(anchor.coor)
    }))
    if(mapLegColor) {
      layout.labels(x = centroids[, 1], y = centroids[, 2], labels = Annot, col = colorspace::darken(colMapping, 0.8), bg = "#eeeeee", r = 0.1, cex = legcex) #colorspace::darken(colMapping, 0.8) #"#eeeeee"
    } else layout.labels(x = centroids[, 1], y = centroids[, 2], labels = Annot, col = legCol, bg = legbg, r = 0.1, cex = legcex) #colorspace::darken(colMapping, 0.8) #"#eeeeee"
    
    #wordcloud::textplot(x = centroids[, 1], y = centroids[, 2], new = F, words = Annot, col = colorspace::darken(Pal, 0.5), bg = "#eeeeee")
  }
  
  if ( (add.legend == TRUE)) {
    xmin <- par("usr")[1]
    xmax <- par("usr")[2]
    ymin <- par("usr")[3]
    ymax <- par("usr")[4]
    
    colMapping <- Col[unique(names(Col))]
    Annot <- names(colMapping)
    labels <- names(Col)
    
    lgd <- legend(x = mean(c(xmin,xmax)), y =  mean(c(ymin,ymax)), legend = Annot, fill = colMapping, cex = 0.5, bty = "n", plot = F)
    
    par(xpd = T, mai = c(0, 0, 0, lgd$rect$w))
    
    legend(x = xmax, y = ymin+lgd$rect$h, legend = Annot, fill = colMapping, cex = 0.5, bty = "n", plot = T)
  }
  
}
######################################################
Plot.cluster.stats.annot <- function(xsce, cell.type.colors=NULL, cell.type="cell.type", clusterLab="cluster", ...) {
  annoMapMatrix <- Extract.cell.annotation.mappings(xsce, ...)
  
  if(is.null(cell.type.colors)) cell.type.colors <- ACTIONet.color.bank0[1:length(xsce@metadata$ClusterCelltypeMapping[,cell.type])]
  par(mfrow=c(6,1), mar=c(3.5, 4, 1.5, 2.5))
  
  boxplot(log(xsce@colData$TotalCounts)~factor(xsce@colData[,clusterLab], levels = unique(annoMapMatrix[,clusterLab])), las=2, ylab="total counts (log)", xlab="")
  
  boxplot(xsce@colData$MitoNonMitoRation~factor(xsce@colData[,clusterLab], levels = unique(annoMapMatrix[,clusterLab])), las=2, ylab="mitochondrial ratio", xlab="")
  
  boxplot(xsce@colData$nGenesDetected~factor(xsce@colData[,clusterLab], levels = unique(annoMapMatrix[,clusterLab])), las=2, ylab="genes detected", xlab="")
  
  barplot(table(xsce@colData[,clusterLab])[annoMapMatrix[,clusterLab]], las=2, col=cell.type.colors[annoMapMatrix[,cell.type]], ylab="cell count", xlab="")
  barplot(rowSums(Table.to.matrix(table(xsce@colData[,clusterLab], xsce@colData$projid))>10)[annoMapMatrix[,clusterLab]], las=2, ylab="individual count (>10cells)", xlab="")
  
  barplot(Counts.to.fractions(Table.to.matrix(table(xsce@colData$projid, xsce@colData[,clusterLab])))[,annoMapMatrix[,clusterLab]], las=2, col=ACTIONet.color.bank, ylab="individuals", xlab="")
}
######################################################
Plot.top.genes <- function(xsce, Col="tomato", mfrowPar=c(4,4), marPar=c(3,8,3,3)) {
  x <- xsce@metadata$cluster.top.genes
  par(mfrow=mfrowPar, mar=marPar)
  for(i in names(x)) barplot(x[[i]][length(x[[i]]):1], las=2, horiz=T, main=i, col=Col)
}
######################################################
Plot.top.genes.all.subclusters <- function(xsce, Col="tomato", ...) {
  #mfrowPar=c(4,4), marPar=c(3,8,3,3)
  for(i in names(xsce@metadata$subclustering.out)) Plot.top.genes(xsce = xsce@metadata$subclustering.out[[i]], Col = Cell.colors[i], ...)
}
######################################################
Plot.subcluster.correlation <- function(xsce, xColors=NULL) {
  if(is.null(xColors)) return("need cell type colors")
  M.sig <- Extract.subcluster.global.signatures(xsce)
  CtLab <- Str.extract(colnames(M.sig), ".", 1, fixed=T)
  Heatmap(cor(M.sig), rect_gp = gpar(col = "black", name="cor"), bottom_annotation = HeatmapAnnotation(cell.type=CtLab, col=list(cell.type=xColors)), name="cor")
}
######################################################
Plot.subcluster.networks <- function(xsce, celltypeLabs="PsychENCODE.celltypes", CPall = ACTIONet.color.bank, ...) {
  L <- split(xsce@colData, xsce@colData[[celltypeLabs]])
  L <- lapply(L, as.data.frame)
  for(i in names(L)) Plot.coords(L[[i]], x = "sub.x", y="sub.y", Col = Vector.to.colors(L[[i]]$sub.cluster, CPall = CPall), add.text = T, ...)
}
######################################################
Plot.subcluster.phylo <- function(xsce, addColor=F) {
  D <- dist(cor(Extract.subcluster.global.signatures(xsce)))
  HC <- hclust(D)
  Col <- "black"
  
  if(addColor) {
    annoMapMatrix <- Extract.cell.annotation.mappings(xsce, varNames = "sub.cluster")
    Col <- Cell.colors[annoMapMatrix[HC$labels, "cell.type"]]
  }
  plot(ape::as.phylo(HC), type="cladogram", tip.color=Col)
}
######################################################
Openfile <- function(File) {
  system(paste("open", File))
}
######################################################
Add.subclustering.analysis <- function(xsce, groupVar="cell.type.annotation", reReduce=TRUE, BatchCorrect = FALSE) {
  temp <- Split.sce.cols(xsce, groupVar)
  if(reReduce) for(i in 1:length(temp)) {
    temp[[i]]@reducedDims[["S_r"]] <- NULL
    if(BatchCorrect) temp[[i]] <- reduce.and.batch.correct.sce.Harmony(temp[[i]], batch.vec = as.character(temp[[i]]$batch)) else temp[[i]] <- reduce.sce(temp[[i]])
  }
  xsce@metadata$subclustering.out <- lapply(temp, function(i) ACTIONet.pattern.analysis(xsce = i))
  try(xsce <- Add.subclustering.metadata(xsce = xsce, groupVar = "PsychENCODE.celltypes"))
  return(xsce)
}
######################################################
Split.sce.cols <- function(InSCE, InClass) {
  temp.sce.split <- lapply(as.character(sort(unique(InSCE@colData[,InClass]))), function(i) InSCE[,InSCE@colData[,InClass]==i])
  names(temp.sce.split) <- as.character(sort(unique(InSCE@colData[,InClass])))
  return(temp.sce.split)
}
######################################################
Add.subclustering.metadata <- function(xsce, groupVar) {
  x <- do.call(rbind, lapply(xsce@metadata$subclustering.out, Extract.subcluster.attributes.from.sce))
  x <- DataFrame(xsce@colData, x[match(colnames(xsce), rownames(x)),])
  x$sub.cluster <- paste(x[,groupVar], x[,"sub.cluster"], sep=".")
  table(x[,groupVar], x$sub.cluster)
  xsce@colData <- x
  return(xsce)
}
######################################################
Extract.subcluster.attributes.from.sce <- function(xsce, preF="sub", extraCols=NULL) {
  temp <- do.call(cbind, vertex.attributes(xsce@metadata$network)[sapply(vertex.attributes(xsce@metadata$network), is.numeric)])
  rownames(temp) <- vertex.attributes(xsce@metadata$network)[["name"]]
  #temp <- cbind(temp, cluster=xsce@metadata$cluster)
  temp <- DataFrame(temp, cluster=xsce@metadata$cluster)
  colnames(temp) <- paste(preF, colnames(temp), sep = ".")
  return(temp)
}
######################################################
Filter.significant.enrichment.map.values <- function(enrichMat, pCut=0.01, keepRelevant=TRUE) {
  pvals <- Bonferronit(convert.z.score(enrichMat))
  enrichMat[pvals>pCut] <- 0
  if(keepRelevant) enrichMat <- enrichMat[,colSums(enrichMat>0)>0]
  return(enrichMat)
}
######################################################
convert.z.score <- function(z, one.sided=NULL) {
  if(is.null(one.sided)) {
    pval = pnorm(-abs(z));
    pval = 2 * pval
  } else if(one.sided=="-") {
    pval = pnorm(z);
  } else {
    pval = pnorm(-z);
  }
  return(pval);
}
######################################################
Bonferronit <- function(x, RemoveNonSig=TRUE, Tresh=0.01) {
  X <- matrix(p.adjust(as.numeric(x), method = "bonferroni"), nrow(x), ncol(x))
  rownames(X) <- rownames(x)
  colnames(X) <- colnames(x)
  #X <- x*length(as.numeric(x))
  #X[X>1] <- 1
  if(RemoveNonSig) X[X>Tresh] <- 1
  return(X)
}
######################################################
Vector.to.colors <- function(inVec, CPall="npg") {
  inVec <- as.character(inVec)
  out <- ggpubr::get_palette(CPall, length(unique(inVec)))[factor(inVec)]
  names(out) <- inVec
  return(out)
}
######################################################
Extract.cell.annotation.mappings <- function(xsce, varNames = c("cluster", "sub.cluster"), celltypeAnnotVar="PsychENCODE.celltypes") {
  CellAnnotations <- unique(xsce@colData[,c(varNames, celltypeAnnotVar)])
  rownames(CellAnnotations) <- CellAnnotations$sub.cluster
  names(CellAnnotations)[names(CellAnnotations)==celltypeAnnotVar] <- "cell.type"
  CellAnnotations <- data.frame(CellAnnotations[order(CellAnnotations$cell.type),])
  return(CellAnnotations)
}
######################################################
Find.outliers <- function(cellinfo, coordLabs2D=c("x","y"), coordLabs3D=c("X", "Y", "Z"), returnMetrics=FALSE) {
  #xmetrics <- lapply(lapply(coordLabs, function(i) lm(cellinfo[,i] ~ ., data = as.data.frame(cellinfo[,coordLabs[!coordLabs%in%i]]))), car::outlierTest)
  #xmetrics[["TD"]] <- car::outlierTest(lm(y~x, data = cellinfo[,c("x","y")]))
  mod <- lm(y ~ x, data = as.data.frame(cellinfo[,coordLabs2D]))
  mod2 <- lm(Y ~ ., data = as.data.frame(cellinfo[,coordLabs3D]))
  temp <- car::outlierTest(mod, n.max=Inf, cutoff=Inf)
  temp2 <- car::outlierTest(mod, n.max=Inf, cutoff=Inf)
  xoutliers <- c(names(which(temp$p<0.01)), names(which(temp2$p<0.01)))
  if(returnMetrics) return(list(mod, mod2))
  return(xoutliers)
}
######################################################
Get.outlier.cells <- function(inCellInfo, cellGrouping="PsychENCODE.celltypes", showPlots=TRUE) {
  temp <- inCellInfo
  temp.cell.info <- split(temp, temp[,cellGrouping])
  Outliers.by.celltype <- lapply(temp.cell.info, Find.outliers)
  
  if(showPlots) {
    for(i in names(temp.cell.info)) {
      Plot.coords(temp.cell.info[[i]], "x", "y", Col = "grey")
      points(as.matrix(temp.cell.info[[i]][Outliers.by.celltype[[i]],c("x", "y")]), col=2, cex=1, pch=20)
    }
  }
  
  toRemove <- unique(unlist(Outliers.by.celltype))
  return(toRemove)
}
######################################################
Get.reLayout.corrdinates <- function(xsce) {
  # get new coordinates using only curated cells
  initial.coordinates = t(scale(xsce@reducedDims[[1]]))
  vis.out = layoutACTIONet(G = get.adjacency(xsce@metadata$network, attr = "weight"), S_r = initial.coordinates, compactness_level = 50, n_epochs = 500)
  xx <- data.frame(do.call(cbind, vis.out[1:2]))
  names(xx) <- paste0("re.", c("x","y", "X", "Y", "Z"))
  gc()
  return(xx)
}
######################################################
Get.column.group.average <- function(x, group, doParallel=FALSE) {
  Ugroups <- sort(names(which(table(group)>1)))
  print(table(group)[Ugroups])
  if(doParallel) out <- do.call(cbind, parallel::mclapply(mc.cores = 8, X=Ugroups, function(i) Matrix::rowMeans(as(x[,group==i], "sparseMatrix")))) else out <- do.call(cbind, lapply(Ugroups, function(i) Matrix::rowMeans(as(x[,group==i], "sparseMatrix"))))
  colnames(out) <- Ugroups
  rownames(out) <- rownames(x)
  return(out)
}
######################################################
Get.column.group.var <- function(x, group, doParallel=FALSE) {
  Ugroups <- sort(names(which(table(group)>1)))
  print(table(group)[Ugroups])
  if(doParallel) out <- do.call(cbind, parallel::mclapply(mc.cores = 8, X=Ugroups, function(i) apply(x[,group==i], 1, var))) else out <- do.call(cbind, lapply(Ugroups, function(i) apply(x[,group==i], 1, var)))
  colnames(out) <- Ugroups
  rownames(out) <- rownames(x)
  return(out)
}
######################################################
Get.column.group.detection.rate <- function(x, group, doParallel=FALSE) {
  Ugroups <- sort(names(which(table(group)>1)))
  print(table(group)[Ugroups])
  if(doParallel) out <- do.call(cbind, parallel::mclapply(mc.cores = 8, X=Ugroups, function(i) Matrix::rowMeans(as(x[,group==i], "sparseMatrix")>0))) else out <- do.call(cbind, lapply(Ugroups, function(i) Matrix::rowMeans(as(x[,group==i], "sparseMatrix")>0)))
  colnames(out) <- Ugroups
  rownames(out) <- rownames(x)
  return(out)
}
######################################################
Get.column.group.count <- function(x, group, doParallel=FALSE) {
  Ugroups <- sort(names(which(table(group)>1)))
  print(table(group)[Ugroups])
  if(doParallel) out <- do.call(cbind, parallel::mclapply(mc.cores = 8, X=Ugroups, function(i) Matrix::rowSums(as(x[,group==i], "sparseMatrix")))) else out <- do.call(cbind, lapply(Ugroups, function(i) Matrix::rowSums(as(x[,group==i], "sparseMatrix"))))
  colnames(out) <- Ugroups
  rownames(out) <- rownames(x)
  return(out)
}
######################################################
Get.pseudobulk.summary <- function(xsce, groupVect, ...) {
  E=Get.column.group.average(xsce@assays[["logcounts"]], group = groupVect, ...)
  C=Get.column.group.count(xsce@assays[["counts"]], group = groupVect, ...)
  V=Get.column.group.var(xsce@assays[["logcounts"]], group = groupVect, ...)
  usage=Get.column.group.detection.rate(xsce@assays[["logcounts"]], group = groupVect, ...)
  PB.sce = SingleCellExperiment(assays = list(E = E, C = C, V = V, usage = usage))
  PB.sce$cell.counts = table(groupVect)[colnames(PB.sce)]
  PB.sce@assays[["expressed"]] <- PB.sce@assays[["usage"]]>0.1
  PB.sce@metadata$expressed.genes <- lapply(data.frame(PB.sce@assays[["expressed"]]), function(i) rownames(PB.sce@assays[["expressed"]])[which(i)])
  return(PB.sce)
}
######################################################
Table.to.matrix <- function(x) {
  return(t(apply(x, 1, function(i) i)))
}
######################################################
Counts.to.fractions <- function(x, bycol=T) {
  if(bycol) return(apply(x, 2, function(i) i/sum(i))) else return(apply(x, 1, function(i) i/sum(i)))
}
######################################################
Aggregate.pairwise.FC <- function(expMat, QueryCol, RefColSet=NULL, sortedFC=F) {
  if(is.null(RefColSet)) RefColSet <- colnames(expMat)[-which(colnames(expMat)==QueryCol)]
  temp <- do.call(cbind, lapply(RefColSet, function(i) expMat[,QueryCol]-expMat[,i]))
  colnames(temp) <- RefColSet
  rownames(temp) <- rownames(expMat)
  FC <- rowSums(temp)
  if(sortedFC) FC <- sort(FC, decreasing = T)
  return(FC)
}
######################################################
Aggregate.pairwise.FC.colPairwise <- function(xM) {
  out <- do.call(cbind, lapply(colnames(xM), function(i) Aggregate.pairwise.FC(xM, QueryCol = i, sortedFC = F)))
  colnames(out) <- colnames(xM)
  return(out)
}
######################################################
Append.to.list <- function(Eli, Li) {
  for(i in 1:length(Eli)) {
    Li[[length(Li)+1]] <- Eli[[i]]
    names(Li)[length(Li)] <- names(Eli)[i]
  }
  return(Li)
}
######################################################
Remove.zero.pvalues <- function(Qmat)  {
  Qmat[Qmat==0] <- min(Qmat[Qmat!=0])*0.1
  return(Qmat)
}
######################################################

sum_counts = function(counts, label, cell_labels){
    # Sums cell-level counts by factors in label vector
    #
    # Args:
    #   counts: a sparse matrix with read counts (gene x cell)
    #   label: variable of interest by which to sum counts
    #   cell_labels: vector of cell labels
    #
    # Returns:
    #   summed counts
    #   number of cells used per summation

    colnames(label) = 'ID'
    label$index = 1
    label$celltype = cell_labels
    label = as.data.frame(pivot_wider(label, values_from = index, names_from = ID))
    label[is.na(label)] = 0
    label$celltype = NULL
    summed_counts = counts%*%as.matrix(label)

    ncells = colSums(label)

    return(list('summed_counts' = summed_counts, 'ncells' = ncells))
}
             
######################################################
normalize.default = function (ace, log_scale = TRUE){ # from ACTIONet V.2.0
    S = SummarizedExperiment::assays(ace)[["counts"]]
    B = rescale.matrix(S, log_scale)
    rownames(B) = rownames(ace)
    colnames(B) = colnames(ace)
    SummarizedExperiment::assays(ace)[["logcounts"]] = B
    return(ace)
}
######################################################
rescale.matrix = function (S, log_scale = FALSE){  # from ACTIONet V.2.0
    if (is.matrix(S)) {
        cs = Matrix::colSums(S)
        cs[cs == 0] = 1
        B = median(cs) * scale(S, center = F, scale = cs)
        if (log_scale == TRUE) {
            B = log1p(B)
        }
    }
    else {
        A = as(S, "dgTMatrix")
        cs = Matrix::colSums(A)
        cs[cs == 0] = 1
        x = median(cs) * (A@x/cs[A@j + 1])
        if (log_scale == TRUE) {
            x = log1p(x)
        }
        B = Matrix::sparseMatrix(i = A@i + 1, j = A@j + 1, x = x, 
            dims = dim(A))
    }
    return(B)
}
                               
######################################################
Extract.tops.by.column <- function(inMatrix, nTop=10, dec=T) {
  out <- lapply(1:ncol(inMatrix), function(i) sort(inMatrix[,i], decreasing = dec)[1:nTop])
  names(out) <- colnames(inMatrix)
  return(out)
}

######################################################
ACTIONet.pattern.analysis <- function(xsce, Renorm=FALSE, doCluster = TRUE, ClustLabel="cluster", ...) {
  if(Renorm) {
    assays(xsce)[[2]] <- NULL
    reducedDim(xsce, type = "S_r") <- NULL
    #reducedDim(xsce)[[1]] <- NULL
    xsce <- reduce.sce(xsce)
  }
  out <- run.ACTIONet(xsce, ...)
  if(doCluster) out <- cluster.ACTIONet(out, annotation.name = ClustLabel)
  xsce@metadata$C <- as(out$unification.out$C.core, "sparseMatrix")
  xsce@metadata$H <- as(out$unification.out$H.core, "sparseMatrix")
  xsce@metadata$W <- out$unification.out$W.core

  xsce@metadata$state.profiles <- out$unification.out$cellstates.core
  colnames(xsce@metadata$state.profiles) <- paste0("P", 1:ncol(xsce@metadata$state.profiles))
  xsce@metadata$state.signatures <- Aggregate.pairwise.FC.colPairwise(xM = xsce@metadata$state.profiles)
  xsce@metadata$signature.top.genes <- Extract.tops.by.column(xsce@metadata$state.signatures)
  xsce@metadata$network <- out$ACTIONet
  V(xsce@metadata$network)$name <- colnames(xsce)
  xsce@metadata$vis.out <- out$vis.out
  xsce@metadata[[ClustLabel]] <- paste0("C", as.character(out$annotations[[ClustLabel]]$Labels))

  xsce@metadata[[paste0(ClustLabel,".profiles")]] <- Get.column.group.average(x = assays(xsce)[["logcounts"]], group = xsce@metadata[[ClustLabel]])
  xsce@metadata[[paste0(ClustLabel,".signatures")]] <- Aggregate.pairwise.FC.colPairwise(xsce@metadata[[paste0(ClustLabel,".profiles")]])
  xsce@metadata[[paste0(ClustLabel,".top.genes")]] <- Extract.tops.by.column(xsce@metadata[[paste0(ClustLabel,".signatures")]])
  xsce@metadata[[paste0(ClustLabel,".cor")]] <- cor(xsce@metadata[[paste0(ClustLabel,".signatures")]])

  DF <- DataFrame(num=1:ncol(xsce), tag=colnames(xsce),
                  x=xsce@metadata$vis.out$coordinates[,1], y=xsce@metadata$vis.out$coordinates[,2],
                  X=xsce@metadata$vis.out$coordinates_3D[,1], Y=xsce@metadata$vis.out$coordinates_3D[,2], Z=xsce@metadata$vis.out$coordinates_3D[,3],
                  color=rgb(xsce@metadata$vis.out$colors),
                  connectivity=V(xsce@metadata$network)$connectivity/max(V(xsce@metadata$network)$connectivity),
                  assignments.dominant=paste0("P",out$unification.out$assignments.core),
                  assignments.dominant.confidence=out$unification.out$assignments.confidence.core,
                  cluster=xsce@metadata$cluster
  )
  xsce@colData <- DataFrame(xsce@colData, DF)
  rm(out)
  gc()
  return(xsce)
}
######################################################
Extract.subcluster.global.signatures <- function(xsce) {
  L <- lapply(names(xsce@metadata$subclustering.out), function(i) Extract.subcluster.sce(xsce, i)@metadata$cluster.profile)
  names(L) <- names(xsce@metadata$subclustering.out)
  CtLab <- rep(names(L), sapply(L, ncol))
  M <- do.call(cbind, L)
  colnames(M) <- paste(CtLab, colnames(M), sep=".")
  M.sig <- Aggregate.pairwise.FC.colPairwise(M)
  return(M.sig)
}
######################################################
Get.top.matching.PanglaoDB <- function(xsignMat) {
  temp <- lapply(PanglaoDB, function(i) Permutation.enrichment.analysis(t(xsignMat), marker.genes = i)$enrichment)
  temp <- Extract.tops.by.column(t(do.call(cbind, temp)))
  return(temp)
}
######################################################
Permutation.enrichment.analysis <- function(x, marker.genes, rand.sample.no = 1000) {
  #require(ACTIONet)
  require(igraph)
  require(Matrix)
  require(stringr)

  if(is.matrix(marker.genes) | is.sparseMatrix(marker.genes)) {
    marker.genes = apply(marker.genes, 2, function(x) rownames(marker.genes)[x > 0])
  }

  archetype.panel <- x

  GS.names = names(marker.genes)
  if (is.null(GS.names)) {
    GS.names = sapply(1:length(GS.names), function(i) sprintf("Celltype %s", i))
  }

  markers.table = do.call(rbind, lapply(names(marker.genes), function(celltype) {
    genes = marker.genes[[celltype]]
    if (length(genes) == 0)
      return(data.frame())


    signed.count = sum(sapply(genes, function(gene) grepl("\\+$|-$", gene)))
    is.signed = signed.count > 0

    if (!is.signed) {
      df = data.frame(Gene = (genes), Direction = +1, Celltype = celltype)
    } else {

      pos.genes = (as.character(sapply(genes[grepl("+", genes, fixed = TRUE)], function(gene) stringr::str_replace(gene,
                                                                                                                   stringr::fixed("+"), ""))))
      neg.genes = (as.character(sapply(genes[grepl("-", genes, fixed = TRUE)], function(gene) stringr::str_replace(gene,
                                                                                                                   stringr::fixed("-"), ""))))

      df = data.frame(Gene = c(pos.genes, neg.genes), Direction = c(rep(+1, length(pos.genes)), rep(-1, length(neg.genes))),
                      Celltype = celltype)
    }
  }))
  markers.table = markers.table[markers.table$Gene %in% colnames(archetype.panel), ]

  if (dim(markers.table)[1] == 0) {
    print("No markers are left")
    return()
  }

  IDX = split(1:dim(markers.table)[1], markers.table$Celltype)

  print("Computing significance scores")
  set.seed(0)
  Z = sapply(IDX, function(idx) {
    markers = (as.character(markers.table$Gene[idx]))
    directions = markers.table$Direction[idx]
    mask = markers %in% colnames(archetype.panel)

    A = as.matrix(archetype.panel[, markers[mask]])
    sgn = as.numeric(directions[mask])
    stat = A %*% sgn

    rand.stats = sapply(1:rand.sample.no, function(i) {
      rand.samples = sample.int(dim(archetype.panel)[2], sum(mask))
      rand.A = as.matrix(archetype.panel[, rand.samples])
      rand.stat = rand.A %*% sgn
    })

    cell.zscores = as.numeric((stat - apply(rand.stats, 1, mean))/apply(rand.stats, 1, sd))

    return(cell.zscores)
  })

  Z[is.na(Z)] = 0
  Labels = colnames(Z)[apply(Z, 1, which.max)]

  #L = names(marker.genes)
  #L.levels = L[L %in% Labels]
  #Labels = match(L, L.levels)
  #names(Labels) = L.levels
  #Labels = factor(Labels, levels = L)
  Labels.conf = apply(Z, 1, max)

  names(Labels) = rownames(archetype.panel)
  names(Labels.conf) = rownames(archetype.panel)
  rownames(Z) = rownames(archetype.panel)

  out.list = list(annotation=data.frame(labels = Labels, labels.confidence = Labels.conf, stringsAsFactors = FALSE), enrichment = Z)

  return(out.list)
}
######################################################

Extract.subcluster.sce <- function(xsce, subLabel){
    x = xsce@metadata$subclustering.out[[subLabel]]
    return(x)
}


Str.extract = function(names){
    x = strsplit(names, '[.]')
    return(unlist(lapply(x, function(i) i[[1]][1])))
}

