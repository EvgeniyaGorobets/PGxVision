% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{Biomarkers}
\alias{Biomarkers}
\title{Drug sensitivity data for genes in the GRCh38.p13 genome assembly}
\format{
A data frame with 7805 rows and 12 variables:
\describe{
  \item{gene}{ENSEMBL gene name}
  \item{compound}{name of compound used in experiment}
  \item{tissue}{name of tissue used in experiment}
  \item{estimate}{The drug sensitivity estimate}
  \item{n}{The number of trials}
  \item{fdr}{false discovery rate (adjusted p-value)}
  \item{mDataType}{the molecular data type, i.e., rna, cnv, etc.}
  \item{gene_seq_start}{the relative starting position of the gene within
  the chromosome}
  \item{gene_seq_end}{the relative ending position of the gene within
  the chromosome}
  \item{strand}{which strand of DNA the gene is located on (+/-)}
  \item{chr}{the chromosome on which the gene is located}
  \item{pvalue}{the p-value of the drug response}
}
}
\source{
Feizi, N., Nair, S. K., Smirnov, P., Beri, G., Eeles, C., Esfahani,
P. N., … Haibe-Kains, B. (2021). PharmacoDB 2.0: Improving scalability and
transparency of in vitro pharmacogenomics analysis. _bioRxiv._
doi:10.1101/2021.09.21.461211
}
\usage{
Biomarkers
}
\description{
A dataset containing drug sensitivity signatures for cancer cell-line
experiments. Each row contains the computed drug sensitivity signature
(estimate, p-value, fdr, n) of a single gene in a single experiment (an
experiment is defined by a tissue, compound/drug, and molecular data type).
}
\keyword{datasets}
