% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GsDenPlot.R
\name{GsDenPlot}
\alias{GsDenPlot}
\title{Plot the density for features}
\usage{
GsDenPlot(
  seu,
  features,
  slot = "data",
  assay = "RNA",
  bias = 1,
  reduction = NULL,
  size = 0.5,
  shape = 16,
  cells = 25,
  pal = "turbo",
  raster = "T",
  ...
)
}
\arguments{
\item{seu}{seurat obj}

\item{slot}{which slot}

\item{assay}{which assay}

\item{bias}{weight of smooth}

\item{reduction}{which  reduction used for knn smoothing}

\item{size}{dot size for plot}

\item{shape}{shape of dot}

\item{cells}{cell number for knn}

\item{pal}{colors}

\item{raster}{logical: raster dot}

\item{feature}{feature for display}
}
\value{
density matirx and gglot object
}
\description{
Plot the density for features
}
\author{
Qiong Zhang
}
