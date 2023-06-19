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
#' @param y vector with 2 values with up and down boundaries
#' @name scaleXY
#' @return rescale vector with value between 0 and 1
#' @author Qiong Zhang
#'
scaleXY <- function(x, y) {
  (x - y[1]) / (y[2] - y[1])
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

#' Transform a list comprised of character vectors with equal lengths to matrix
#'
#' @param vectorList Vector list
#' @name sameSizeVectorList2Matrix
#' @return Matrix contains character
#' @author Qiong Zhang
#'
sameSizeVectorList2Matrix <- function(vectorList){
  sm_i<-NULL
  sm_j<-NULL
  sm_x<-NULL
  for (k in 1:length(vectorList)) {
    sm_i <- c(sm_i,rep(k,length(vectorList[[k]]@i)))
    sm_j <- c(sm_j,vectorList[[k]]@i)
    sm_x <- c(sm_x,vectorList[[k]]@x)
  }
  return (sparseMatrix(i=sm_i,j=sm_j,x=sm_x,dims=c(length(vectorList),vectorList[[1]]@length)))
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


#' Calculate normalized enrichment score for gene sets. To overcome the dropout effects in scRNA-seq data, imputate gene ranking via gene co-expression modules detected by affinity propagation (AP).
#'
#' @param expMatRank Ranking matrix of gene (row) X cell (col) matrix
#' @param geneSets Genesets
#' @param tfidfMatrix tfidfMatrix matrix
#' @param imputation logical value for gene ranking imputation
#' @param retriRatio Drop out ratio of genes in given gene set.
#' @param retriCellRatio The ratio of cells containing drop out genes.
#' @name calGSnorm Calculate normalized gene signatures for gene sets
#' @return Signature matrix for gene sets
#' @author Qiong Zhang
#'
calGSnorm <- function(expMatRank, geneSets, tfidfMatrix, imputation = F, scale = scale, retriRatio = 0.1, retriCellRatio = 0.1) {
  .gene.all <- rownames(expMatRank)
  .m <- calRange(geneSets, .gene.all)
  .cell.n <- ncol(expMatRank)
  .score.lst <- lapply(seq_along(geneSets), function(x){
    .genes.a <- geneSets[[x]]
    .gene.rev.index <- grepl('-$', .genes.a, perl=TRUE)
    .gene.index <- !.gene.rev.index
    .genes <- .genes.a[.gene.index]
    .gene.n <- length(.genes)
    .gs.ind <- match(.genes, .gene.all)
    .rankMat.sub <- expMatRank[.gs.ind, , drop = FALSE]
    if (isTRUE(imputation)){
      .genes.NonExp <- matrixStats::colCounts(.rankMat.sub, value = 1)
      .genes.NonExp.ratio <- .genes.NonExp / .gene.n
      .retri.index <- which(.genes.NonExp.ratio <= retriRatio)
      if (length(.retri.index) > .cell.n * retriCellRatio ){
        .matrix.subset <- tfidfMatrix[.gs.ind, seq(1, .cell.n, 10)]
        .gene.sim <- proxyC::simil(.matrix.subset, .matrix.subset, method = "cosine")
        .gene.sim.apcluster <- apcluster::apcluster(as.matrix(.gene.sim))
        .gene.sim.apcluster.info <- lapply(.gene.sim.apcluster@clusters, function(x) if(length(x)>4) names(x))
        .gene.sim.apcluster.info.filter <- .gene.sim.apcluster.info[-which(sapply(.gene.sim.apcluster.info, is.null))]
        if (length(.gene.sim.apcluster.info.filter) > 0 ){
          .genes.NonExp.all <- matrixStats::colCounts(expMatRank, value = 1)
          .gene.sim.rep.index <- match(unlist(.gene.sim.apcluster.info.filter), .genes)
          .retrieval.matrix <- do.call(rbind, lapply(seq_along(.gene.sim.apcluster.info.filter), function(x){
            .cluster.genes.tmp.name <- .gene.sim.apcluster.info.filter[[x]]
            .cluster.genes.tmp.name.len <- length(.cluster.genes.tmp.name)
            .cluster.genes.tmp.index <- match(.cluster.genes.tmp.name, .genes)
            .cluster.genes.tmp.matrix <- .rankMat.sub[.cluster.genes.tmp.index, .retri.index]
            #####
            .genes.NonExp.cluster <- matrixStats::colCounts(.cluster.genes.tmp.matrix, value = 1)
            .retri.index.cluster <- which(.genes.NonExp.cluster / .cluster.genes.tmp.name.len <= retriRatio)
            if (length(.retri.index.cluster) > 0 ){
              .cluster.genes.tmp.matrix.sub <- .cluster.genes.tmp.matrix[, .retri.index.cluster]
              .rank.matrix.rep.p <- .cluster.genes.tmp.matrix.sub == 1
              .rank.matrix.rep.v <-  t(replicate(.cluster.genes.tmp.name.len, .genes.NonExp.all[.retri.index.cluster]))
              .rank.matrix.rep.new <- .cluster.genes.tmp.matrix.sub + .rank.matrix.rep.p * .rank.matrix.rep.v
              .cluster.genes.tmp.matrix[, .retri.index.cluster] <- .rank.matrix.rep.new
            }
            return(.cluster.genes.tmp.matrix)
          }))
          .rankMat.sub[.gene.sim.rep.index, .retri.index] <- .retrieval.matrix
        } else {
          warning("Failure when retrieving gene ranking due to the lack of gene module")
        }
      }
    }
    .gs.rank <- matrixStats::colMeans2(.rankMat.sub)
    .score <- scaleXY(.gs.rank, .m[,x])
    if (any(.gene.rev.index)){
      .genes.rev <- .genes.a[.gene.rev.index]
      .gs.ind.rev <- match(.genes.rev, .gene.all)
      .rankMat.sub.rev <- expMatRank[.gs.ind.rev, , drop = FALSE]
      .gs.rank.rev <- matrixStats::colMeans2(.rankMat.sub.rev)
      .score.rev <- scaleXY(.gs.rank, .m[,x])/2
      if (isTRUE(scale))  .score <- scale01(.score - .score.rev)
    }
    return(.score)
  })
  .score.lst
}


#' Calculate normalized enrichment score for gene sets. To overcome the dropout effects in scRNA-seq data, imputate gene ranking via gene co-expression modules detected by affinity propagation (AP).
#'
#' @param expMatRank Ranking matrix of gene (row) X cell (col) matrix
#' @param geneSets Genesets
#' @param tfidfMatrix tfidfMatrix matrix
#' @param imputation logical value for gene ranking imputation
#' @param retriRatio Drop out ratio of genes in given gene set.
#' @param retriCellRatio The ratio of cells containing drop out genes.
#' @name calGSEscore Calculate normalized gene signatures for gene sets
#' @return Signatures gene sets
#' @author Qiong Zhang
#'
calGSEscore <- function (expMatRank, geneSets, tfidfMatrix, imputation = F, scale = F, retriRatio = 0.1, retriCellRatio = 0.1) {
  .GSscore.neg <- 0
  .geneAll <- rownames(expMatRank)
  .geneNum.all <- length(.geneAll)
  .GSscore.tot <- calGSnorm(expMatRank, geneSets, tfidfMatrix, imputation = imputation, scale = scale, retriRatio = retriRatio, retriCellRatio = retriCellRatio)
  return(list(GSig = .GSscore.tot))
}

#' Calculate normalized enrichment score for gene sets. To overcome the dropout effects in scRNA-seq data, imputate gene ranking via gene co-expression modules detected by affinity propagation (AP).
#'
#' @param Seuobj Seurat object
#' @param assay Which assay used in the analysis
#' @param slot Which slot used in the analysis
#' @param geneSets Gene sets, genes with negative effects should be label with "-" suffix, e.g. c("CD8A", "CD8B", "CD4-", "KLRF1-")
#' @param filterRM Filter the ribosome and mito genes before analysis
#' @param normalByRow normalize the importance of genes across cells
#' @param imputation logical value for gene ranking imputation
#' @param retriRatio Drop out ratio of genes in given gene set.
#' @param retriCellRatio The ratio of cells containing drop out genes.
#' @name GSignatures Calculate gene set signatures in scRNA-seq data
#' @return Signatures for gene sets
#' @export
#' @author Qiong Zhang
#'
GSignatures <- function(Seuobj, assay = "RNA", slot = "counts", geneSets, scale = F, filterRM = F, normalByRow = F, imputation = F, retriRatio = 0, retriCellRatio = 0.1){
  require.pkgs <- c("tidyverse", "Seurat", "textmineR", "proxyC", "apcluster")
  invisible(lapply(require.pkgs, function(x) require(x, character.only = T, quietly = T)))
  .expMat <- assay2mtx(Seuobj, assay = assay, slot = slot)
  .genes.all <- rownames(.expMat)
  if(isTRUE(filterRM)){
    .qual_qc_terms <- tibble::tibble(word_reg = c("^RP[SL]", "^MT-", "^HB[^(P)]"), term_name = c("percent_ribo", "percent_mito", "percent_hb"))
    .genes.filter.preset <- grep(str_c(.qual_qc_terms$word_reg, collapse = "|"), .genes.all, value = T)
    .expMat <- .expMat[!(.genes.all %in% .genes.filter.preset), ]
    .expMat <- .expMat[ , matrixStats::colSums2(.expMat) > 100]
    .genes.all <- rownames(.expMat)
  }
  .geneSets.n <- length(geneSets)
  .gs.name <- names(geneSets)
  .GSES.name <- .gs.name
  .geneSet.diff <- setdiff(gsub("-$", "", unlist(geneSets)), .genes.all)
  if (length(.geneSet.diff) > 0){
    .geneSet.diff.string <- paste(.geneSet.diff, collapse = ', ')
    .geneSet.diff.string.err <- paste("Genes missing: ", .geneSet.diff.string)
    stop( .geneSet.diff.string.err )
  }
  .cell.all <- colnames(.expMat)
  .cell.n <- ncol(.expMat)
  .gene.n <- nrow(.expMat)
  .expMat.t <- t(.expMat)
  .expMat.t.freq.mat <- textmineR::TermDocFreq(dtm = .expMat.t)
  .expMat.t.freq.mat$idf[.expMat.t.freq.mat$idf == 0 ] <- 1e-10
  .expMat.t.freq.mat$idf[ is.infinite(.expMat.t.freq.mat$idf) ] <- 0
  .expMat.tf_idf.NormByCell <- t(.expMat.t / Matrix::rowSums(.expMat.t)) * .expMat.t.freq.mat$idf
  if (isTRUE(normalByRow)){
    .expMat.tf_idf.NormByCell.no0 <- .expMat.tf_idf.NormByCell > 0
    .expMat.tf_idf.NormByCell.sum <- Matrix::rowSums(.expMat.tf_idf.NormByCell)
    .expMat.tf_idf.NormByCell[.expMat.tf_idf.NormByCell.no0] <- .expMat.tf_idf.NormByCell[.expMat.tf_idf.NormByCell.no0]/(replicate(.cell.n, .expMat.tf_idf.NormByCell.sum)[as.matrix(.expMat.tf_idf.NormByCell.no0)])
  }
  .expMat.tf_idf.NormByCell.rank <- matrixStats::colRanks(as.matrix(.expMat.tf_idf.NormByCell),
                                                          ties.method = "min",
                                                          preserveShape = TRUE)

  .GSignatures <- calGSEscore(.expMat.tf_idf.NormByCell.rank, geneSets, .expMat.tf_idf.NormByCell, scale = scale, imputation = imputation, retriRatio = retriRatio, retriCellRatio = retriCellRatio)
  names(.GSignatures$GSig) <- .GSES.name
  return(.GSES)
}



