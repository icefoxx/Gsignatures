% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Gsignatures.R
\name{calGSEscore Calculate normalized gene signatures for gene sets}
\alias{calGSEscore Calculate normalized gene signatures for gene sets}
\alias{calGSEscore}
\title{Calculate normalized enrichment score for gene sets. To overcome the dropout effects in scRNA-seq data, imputate gene ranking via gene co-expression modules detected by affinity propagation (AP).}
\usage{
calGSEscore(
  expMatRank,
  geneSets,
  tfidfMatrix,
  imputation = F,
  scale = F,
  retriRatio = 0.1,
  retriCellRatio = 0.1
)
}
\arguments{
\item{expMatRank}{Ranking matrix of gene (row) X cell (col) matrix}

\item{geneSets}{Genesets}

\item{tfidfMatrix}{tfidfMatrix matrix}

\item{imputation}{logical value for gene ranking imputation}

\item{retriRatio}{Drop out ratio of genes in given gene set.}

\item{retriCellRatio}{The ratio of cells containing drop out genes.}
}
\value{
Signatures gene sets
}
\description{
Calculate normalized enrichment score for gene sets. To overcome the dropout effects in scRNA-seq data, imputate gene ranking via gene co-expression modules detected by affinity propagation (AP).
}
\author{
Qiong Zhang
}
