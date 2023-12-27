
#' Print message
#'
#' @param ... Character string to print
#' @param time Logical to print time
#' @param verbose Logical to print message
#' @name print_message
#' @author Qiong Zhang
print_message <- function(..., time = T, verbose = T){
  if (verbose){
    if (time){
      message(Sys.time() , ": ",  ...)
    } else {
      message(...)
    }
  }
}


#' Fetch signal genes existed in expression matrix
#'
#'
#' @name getExistGenes
#' @author Qiong Zhang
#' @param expM Seurat object or expression matrix
#' @param geneSets which feature(s)
#' @return GeneSets contain genes in the expression matrix
#' @export
#'
getExistGenes <- function(expm, geneSets, trim = T, ratio_warn = T, topoint = T){

  .genes.all <- rownames(expm)
  .genes.raw <- unlist(geneSets)
  .genes.inq <- .genes.raw
  if (any(stringr::str_detect(.genes.raw, '-$'))) .genes.inq <- gsub("-$", "", .genes.raw)
  if (isTRUE(trim)) .genes.inq <- gsub('\\.\\d+', '', .genes.inq)
  .genes.inq.f <- .genes.inq %in% .genes.all
  if (!all(.genes.inq.f)) {
    if (isTRUE(ratio_warn)) {
      .ratio.g <- sum(.genes.inq.f) / length(.genes.inq.f) * 100
      if (.ratio.g < 50) warning(stringr::str_glue('{.ratio.g}% genes in the geneSets do NOT exist in the expression matrix'))
    }
  }
  .geneSets.loci.end <- cumsum(sapply(geneSets, length))
  .geneSets.loci.start <- c(1, .geneSets.loci.end[-length(.geneSets.loci.end)]+1)
  .geneSets.f <- lapply(seq_along(.geneSets.loci.end), function(x) {
    .gene.region <- c(.geneSets.loci.start[x] : .geneSets.loci.end[x])
    .genes <- .genes.raw[.gene.region][.genes.inq.f[.gene.region]]
    })
  if (isTRUE(topoint)) {
    names(geneSets) <- stringr::str_replace_all(names(geneSets), pattern = '_', '.')
    message("Feature names with underscores ('_') have been replaced with dot ('.')")
  }
  names(.geneSets.f) <- names(geneSets)
  .geneSets.f.f <- Filter(length, .geneSets.f)
  if (isTRUE(ratio_warn)) {
    .geneSets.missing <- setdiff(names(geneSets), names(.geneSets.f.f))
    if (length(.geneSets.missing) > 0 ) warning(stringr::str_glue('Genes of {.geneSets.missing} do NOT exist in the expression matrix'))
  }
  return(.geneSets.f.f)
}


#' Extract expression matrix from Seurat object or matrix
#'
#' This function extracts the expression matrix from either a Seurat object or a matrix object. If a Seurat object is provided, the function can optionally specify the assay and slot to be used for extracting the expression matrix. If a matrix is provided, the function will create a Seurat object with the provided matrix and default metadata fields. The function can also specify the minimum number of cells and features required to keep genes in the expression matrix.
#'
#' @name assay2mtx
#' @author Qiong Zhang
#' @param object A Seurat object or matrix containing expression data
#' @param assay A character string specifying the assay name to extract expression data from the Seurat object
#' @param slot A character string specifying the slot name to extract expression data from the Seurat object
#' @param min.cells An integer specifying the minimum number of cells required to keep genes in the expression matrix
#' @param min.features An integer specifying the minimum number of features required to keep genes in the expression matrix
#' @return A matrix containing the expression data, with gene names as row names and cell barcodes as column names
#' @examples
#' expmtx <- assay2mtx(object = seurat_object, assay = "RNA", slot = "counts", min.cells = 10, min.features = 5)
#' expmtx <- assay2mtx(object = expression_matrix, min.cells = 10, min.features = 5)
#'
assay2mtx <- function(object, assay = NULL, slot = "counts", min.cells = 3, min.features = 1, update = F) {
  if (any(methods::is(object) %in% "Seurat")) {
    assay <- if (is.null(assay)) Seurat::DefaultAssay(object) else assay
    object <- if (isTRUE(update)) Seurat::UpdateSeuratObject(object) else object
  } else {
    object <- Seurat::CreateSeuratObject(counts = object, project = "GSES", assay = assay, min.cells = min.cells, min.features = min.features)
    assay <- "RNA"
  }
  Seurat::GetAssayData(object, assay = assay, slot = slot)
}

#' rescale a vector between 0 (lowest) and 1 (highest)
#'
#' @param x vector
#' @param y vector with 2 values with max and min value
#' @name scale01
#' @return rescale vector with value between 0 and 1
#' @author Qiong Zhang
#'
scale01 <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

#' rescale a vector between 0 (lowest) and 1 (highest)
#'
#' @param x vector
#' @param y vector with 2 values with up and down boundaries
#' @name scaleXY
#' @return rescale vector with value between 0 and 1
#' @author Qiong Zhang
#'
scaleXY <- function(x, y) {
  (x - y[1]) / (y[2] - y[1])
}

#' Calculate the up- and down- range of ranking for gene sets with gene background
#'
#' @param geneSets Genesets
#' @param geneAll Background genes
#' @name calRange
#' @return Up- and down- range of ranking for gene sets
#' @author Qiong Zhang
#'
calRange <- function(geneSets, geneAll) {
  .geneSets.name <- names(geneSets)
  .ranges <- do.call(cbind, lapply(geneSets, function(x){
    .gene.num.gs <- length(x)
    .gene.num.al <- length(geneAll)
    .z <- ceiling((.gene.num.gs + 1)/2)
    return(c(.z,  .gene.num.al - .z + 1 ))
  }))
  .ranges
}


DateMes <- function(mes){
  message(paste0("[", Sys.time(), "] ", mes))
}


getDatatype <- function(x){
  if(inherits(x, 'Seurat')) return('seu')
  if(inherits(x, 'data.frame')) return('df')
  if (inherits(x, 'matrix')) return('mtx')
  if (inherits(x, 'character')) return('cha')
}

inRange <- function(x,y){
  x >= y[1] & x <= y[2]
}


smKNNCore <- function(s, n, v, w){
    .sn <- c(s,n)
    .w <- (1-w)^(0:length(n))
    sum(.w * v[.sn])/sum(.w)
}


smoothKNN <- function(expM, nn, wt = 0.1, bidirect = F) {
  .knn.scores <- vapply(ncol(expM), FUN.VALUE=numeric(nrow(expM)), FUN=function(f) {
    .v <- expM[,f]
    .knn.s <- vapply(X = seq_len(nrow(nn)),  FUN.VALUE = numeric(1), FUN = function(x) {smKNN(x, nn[x,], .v, wt) })
    if (isFALSE(bidirect)) .knn.s <- pmax(.knn.s, .v)
    return(.knn.s)
  })
  return(.knn.scores)
}


smoothScores <- function(seu, features, assay = 'RNA', slot = "data",  k = 10,
                         wt = 0.1, bidirect = FALSE, reduction = 'umap', dim = 2,
                         method = 'KNN') {
  if (!inRange(wt, c(0,1))) stop("wt MUST be between 0 and 1")
  .k <- ifelse(k<=0, 1, k)
  .data.type.seu <- getDatatype(seu)
  .data.type.redct <- getDatatype(reduction)
  if(.data.type.seu == 'seu') {
    DefaultAssay(seu) <- assay
    .expM <- Seurat::FetchData(seu, vars = features, slot = slot)
    .redct <-Seurat::Reductions(seu)
    if (.data.type.redct == 'cha'){
      if (!reduction %in% .redct) stop(DateMes(stringr::str_glue("{reduction} missing in this object")))
      .xy <- Seurat::Embeddings(seu, reduction = reduction)[, 1:dim]
    } else if (.data.type.redct %in% c('mtx', 'df')) {
      .xy <- reduction[, 1:dim]
    } else {
      stop(DateMes('reduction MUST be a matrix (data.frame) contains cell coordinate information format, or reduction in seurat object'))
    }
  } else if (.data.type.seu %in% c('mtx', 'df')){
    .expM <- seu
    .xy <- reduction[, 1:dim]
  } else {
    stop(DateMes('seu MUST be seurat object or expression matrix with data.frame format'))
  }
  .cell.n <- nrow(.expM)
  if (.cell.n <= .k) .k <- .cell.n - 2
  if (.cell.n > 3) {
    .nn <- BiocNeighbors::findKNN(.xy, k = .k)
    .smooth.scores <- get(paste0('smooth', method))(.expM, .nn$index, wt = wt, bidirect = bidirect)
    .smooth.scores <- as.data.frame(.smooth.scores)
    colnames(.smooth.scores) <- colnames(.expM)
    rownames(.smooth.scores) <- rownames(.expM)
  } else {
    .smooth.scores <- .expM
  }
  return(.smooth.scores)
}


calPurity <- function(expS, gs = 0.6, shift = 0.2, ratio = 2, flag = c('other', 'unclear', 'pure'), neg = F, tofactor = F) {
  .gs <- gs
  if (isTRUE(neg)) .gs <- gs - ratio * shift
  .neg <- shift
  .neg <- ifelse(.neg < 0, 0, .neg)
  .pur <- rep(flag[2], length(expS))
  .pur[which(expS > .gs)] <- flag[3]
  .pur[which(expS < .neg)] <- flag[1]
  if (isTRUE(tofactor)) .pur <- factor(.pur, levels = flag)
  return(.pur)
}




