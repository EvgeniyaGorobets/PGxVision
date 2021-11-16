#' Paxlitaxel response in BRCA PDXs
#'
#' A dataset containing the response of different BRCA tumours to paxlitaxel.
#' The experiment tumours are patient-derived xenographs.
#'
#' @format A data frame with 38 rows and 3 variables:
#' \describe{
#'   \item{tumour}{The name of the PDX tumour}
#'   \item{ODC1}{???}
#'   \item{angle}{The angle between the treatment and control}
#' }
#' @source Xeva. R package version 1.10.0.
"BRCA.PDXE.paxlitaxel.response"


#' Drug sensitivity data for genes in the GRCh38.p13 genome assembly
#'
#' A dataset containing the prices and other attributes of almost 54,000
#' diamonds.
#'
#' @format A data frame with 7805 rows and 12 variables:
#' \describe{
#'   \item{gene}{ENSEMBL gene name}
#'   \item{compound}{name of compound used in experiment}
#'   \item{tissue}{name of tissue used in experiment}
#'   \item{estimate}{???}
#'   \item{n}{???}
#'   \item{fdr}{false discovery rate (adjusted p-value)}
#'   \item{mDataType}{the molecular data type, i.e., rna, etc.}
#'   \item{gene_seq_start}{the relative starting position of the gene within
#'   the chromosome}
#'   \item{gene_seq_end}{the relative ending position of the gene within
#'   the chromosome}
#'   \item{strand}{which strand of DNA the gene is located on (+/-)}
#'   \item{chr}{the chromosome on which the gene is located}
#'   \item{pvalue}{the p-value of the drug response}
#' }
#' @source \url{https://pharmacodb.pmgenomics.ca/}
"Biomarkers"


#' Information about chromosomes in the GRCh38.p13 genome assembly
#'
#' A dataset containing the prices and other attributes of almost 54,000
#' diamonds.
#'
#' @format A data frame with 25 rows and 6 variables:
#' \describe{
#'   \item{# unit-name}{???}
#'   \item{chrName}{the name of the chromosome (all, 1-22, X, Y)}
#'   \item{molecule-type/loc}{the molecule type ("chromosome" or "all")}
#'   \item{sequence-type}{???}
#'   \item{statistic}{the chromosome statistic measured in this row}
#'   \item{chrLength}{the length of the chromosome (# nucleotides)}
#' }
#' @source \url{https://www.gencodegenes.org/pages/gencode.html}
"GRCh38.p13.Assembly"

#' A set of 10 'GO cellular compartment' gene sets
#'
#' A dataset of all gene sets in the category 'GO cellular compartment' that
#' contain the gene ENSG00000012124. Each row represents a single gene in one
#' of the 10 gene sets.
#'
#' @format A data frame with 8187 rows and 2 variables:
#' \describe{
#'   \item{gs_id}{GO gene set ID}
#'   \item{ensembl_gene}{ENSEMBL gene ID}
#' }
#' @source \url{
#' https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=GO:CC}
"TestGeneSets"

# [END]
