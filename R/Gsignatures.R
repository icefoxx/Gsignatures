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
calGSnorm <- function(expMatRank, geneSets, tfidfMatrix, imputation = F, scale = F, retriRatio = 0.1, retriCellRatio = 0.1, decay = 0.5) {
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
      .genes.rev <- gsub("-$", "", .genes.a[.gene.rev.index])
      .gs.ind.rev <- match(.genes.rev, .gene.all)
      .rankMat.sub.rev <- expMatRank[.gs.ind.rev, , drop = FALSE]
      .gs.rank.rev <- matrixStats::colMeans2(.rankMat.sub.rev)
      .score.rev <- scaleXY(.gs.rank.rev, .m[,x]) * decay
      .score <- .score - .score.rev
    }
    if (isTRUE(scale))  .score <- scale01(.score)
    return(.score)
  })
  .score.lst
}


#' Calculate weighted geneset scoring method
#'
#' @param expMatRank Ranking matrix of gene (row) X cell (col) matrix
#' @param geneSets Genesets
#' @param binwidth binwidth for scoring scale
#' @param decay penalty for genes below the mean ranking
#' @return Signatures of gene sets in list object
#' @author Qiong Zhang
#'
weightedBinScore <- function(expMatRank, geneSets, binwidth = 0.2, decay = 0.5) {
  .cell.n <- ncol(expMatRank)
  .gene.all <- rownames(expMatRank)
  .g.max <- length(.gene.all)
  .g.non <- matrixStats::colCounts(expMatRank, value = 1)
#  .g.mean <- matrixStats::colMeans2(expMatRank)
  .g.mean <- (2 * .g.max - .g.non)/2
  .g.bin.coe <- seq(-1, 1, binwidth)
  .g.bin.n <- length(.g.bin.coe)
  .g.bin.coe.m <- median(.g.bin.coe)
  .penalty.index <- which(.g.bin.coe < .g.bin.coe.m)
  .g.bin.coe[.penalty.index] <- .g.bin.coe[.penalty.index] * decay
  .g.c.bin <- seq(0, .g.max, length.out = .g.bin.n + 1)
  .score.lst <- lapply(seq_along(geneSets), function(x){
    .genes.a <- geneSets[[x]]
    .gene.rev.index <- grepl('-$', .genes.a, perl=TRUE)
    .gene.index <- !.gene.rev.index
    .genes <- .genes.a[.gene.index]
    .gene.n <- length(.genes)
    .gs.ind <- match(.genes, .gene.all)
    .rankMat.sub <- expMatRank[.gs.ind, , drop = FALSE]
    .score <- sapply(1:.cell.n, function(y){
      .g.r <- .rankMat.sub[,y]
      .g.in.rank.bins <- cut(.g.r, .g.c.bin)
      .g.in.rank.bins.l <- levels(.g.in.rank.bins)
      .score.s.p <- .g.bin.coe[match(.g.in.rank.bins, .g.in.rank.bins.l)]
      .score.s.r <- .g.r - .g.mean[y]
      .score.s <- sum(.score.s.p * abs(.score.s.r))/.gene.n
      .score.max <- .g.max - .g.mean[y] - .gene.n/2
      .score.min <- .g.mean[y] - 1
      .score.scale <- .score.max + .score.min
      .score.s <- .score.s / .score.scale
      return(.score.s)
    })
  })
  return(.score.lst)
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
calGSEscore <- function (expMatRank, geneSets, tfidfMatrix, imputation = F, scale = F, retriRatio = 0.1, retriCellRatio = 0.1, decay = 1, method = 'GSnorm') {
  .GSscore.neg <- 0
  .geneAll <- rownames(expMatRank)
  .geneNum.all <- length(.geneAll)
  if (method == 'GSnorm'){
    .GSscore.tot <- calGSnorm(expMatRank, geneSets, tfidfMatrix, imputation = imputation, scale = scale, retriRatio = retriRatio, retriCellRatio = retriCellRatio, decay = decay)
  } else if (method == 'weightedBin'){
    .GSscore.tot <- weightedBinScore(expMatRank, geneSets)
  }
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
#' @param method Which method used for Gsignature calculate: 'GSnorm' or 'weightedBin'
#' @name GSignatures Calculate gene set signatures in scRNA-seq data
#' @return seu obj with new assay named as 'GSignatures'
#' @export
#' @author Qiong Zhang
#'
GSignatures <- function(Seuobj, geneSets, assay = "RNA", slot = "counts", add = F, scale = F, filterRM = F, normalByRow = F, imputation = F, retriRatio = 0, retriCellRatio = 0.1, ignore = T, decay = 0.5, method = 'GSnorm', binwidth = 0.2){
  require.pkgs <- c("tidyverse", "Seurat", "textmineR", "proxyC", "apcluster")
  invisible(lapply(require.pkgs, function(x) require(x, character.only = T, quietly = T)))
  .expMat <- assay2mtx(Seuobj, assay = assay, slot = slot)
  .genes.all <- rownames(.expMat)
  if(isTRUE(filterRM)){
    .qual_qc_terms <- tibble::tibble(word_reg = c("^RP[SL]", "^MT-", "^HB[^(P)]"), term_name = c("percent_ribo", "percent_mito", "percent_hb"))
    .genes.filter.preset <- grep(stringr::str_c(.qual_qc_terms$word_reg, collapse = "|"), .genes.all, value = T)
    .expMat <- .expMat[!(.genes.all %in% .genes.filter.preset), ]
#    .expMat <- .expMat[ , matrixStats::colSums2(.expMat) > 100]
    .genes.all <- rownames(.expMat)
  }
  .geneSets.n <- length(geneSets)
  .gs.name <- gsub("[_.]", "-", names(geneSets))
  .GSES.name <- .gs.name
  .geneSet.diff <- setdiff(gsub("-$", "", unlist(geneSets)), .genes.all)
  if (length(.geneSet.diff) > 0 & isFALSE(ignore)){
    .geneSet.diff.string <- paste(.geneSet.diff, collapse = ', ')
    .geneSet.diff.string.err <- paste("Genes missing: ", .geneSet.diff.string)
    stop( .geneSet.diff.string.err )
  }
  .genesets <- getExistGenes(.expMat, geneSets = geneSets)
  .cell.all <- colnames(.expMat)
  .cell.n <- ncol(.expMat)
  .gene.n <- nrow(.expMat)
  .expMat.t <- Matrix::t(.expMat)
  .expMat.t.freq.mat <- textmineR::TermDocFreq(dtm = .expMat.t)
  .expMat.t.freq.mat$idf[.expMat.t.freq.mat$idf == 0 ] <- 1e-10
  .expMat.t.freq.mat$idf[ is.infinite(.expMat.t.freq.mat$idf) ] <- 0
  .expMat.tf_idf.NormByCell <- Matrix::t(.expMat.t / Matrix::rowSums(.expMat.t)) * .expMat.t.freq.mat$idf
  if (isTRUE(normalByRow)){
    .expMat.tf_idf.NormByCell.no0 <- .expMat.tf_idf.NormByCell > 0
    .expMat.tf_idf.NormByCell.sum <- Matrix::rowSums(.expMat.tf_idf.NormByCell)
    .expMat.tf_idf.NormByCell[.expMat.tf_idf.NormByCell.no0] <- .expMat.tf_idf.NormByCell[.expMat.tf_idf.NormByCell.no0]/(replicate(.cell.n, .expMat.tf_idf.NormByCell.sum)[as.matrix(.expMat.tf_idf.NormByCell.no0)])
  }
  .expMat.tf_idf.NormByCell.rank <- matrixStats::colRanks(as.matrix(.expMat.tf_idf.NormByCell),
                                                          ties.method = "min",
                                                          preserveShape = TRUE)
  rownames(.expMat.tf_idf.NormByCell.rank) <- rownames(.expMat.tf_idf.NormByCell)
  .GSignatures <- calGSEscore(.expMat.tf_idf.NormByCell.rank, .genesets, .expMat.tf_idf.NormByCell, scale = scale, imputation = imputation, retriRatio = retriRatio, retriCellRatio = retriCellRatio, decay = decay,method = method)
  .Sample.GSignatures <- do.call(rbind, .GSignatures$GSig)
  rownames(.Sample.GSignatures) <- .gs.name
  colnames(.Sample.GSignatures) <- colnames(.expMat)
  .Sample.GSignatures[.Sample.GSignatures < 0] <- 0
  if (isTRUE(add) & "GSignatures" %in% Seurat::Assays(Seuobj)){
    .tmp.gs <- assay2mtx(Seuobj, assay = "GSignatures", slot = "data")
    .tmp.gs.name <- rownames(.tmp.gs)
    .index.drop <- match(.gs.name, .tmp.gs.name)
    .index.drop.na <- is.na(.index.drop)
    .index.drop <- .index.drop[!is.na(.index.drop)]
    if (length(.index.drop) > 0) .tmp.gs <- .tmp.gs[-.index.drop, ]
    .Sample.GSignatures <- rbind(.tmp.gs, .Sample.GSignatures)
  }
  Seuobj[["GSignatures"]] <- SeuratObject::CreateAssayObject(counts = .Sample.GSignatures)
  return(Seuobj)
}

