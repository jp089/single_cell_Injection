% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_query.r
\name{load_query}
\alias{load_query}
\title{Loading query}
\usage{
load_query(query, adt, list_of_interest, data_name)
}
\arguments{
\item{query}{Path to the directorey containing matrix.mtx,
genes.tsv, and barcodes.tsv files provided by 10X.}

\item{adt}{Path to the input file (tab separated format)}

\item{list_of_interest}{Path to the input file (tab separated format)}

\item{data_name}{name of the data query}
}
\value{
A Seurat object, with adt data as an 'assay'
}
\description{
This function loads gene expression data (query) by
using Read10X (Seurat package), and prepeared to inject
on the reference. If ADT data is availible, it could be added.
If it exists a list of interest, this could be also take
into account to subset both gene expression and ADT data.
Finally, a 'data name' is required. Cells in both gene
expression and ADT data need to be in the same order.
}
