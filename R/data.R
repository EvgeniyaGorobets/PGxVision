#' Paxlitaxel response in BRCA PDXs
#'
#' A dataset containing the response of different BRCA tumours to paxlitaxel,
#' as well as the levels of ODC1 transcripts in each tumour. The tumours are
#' patient-derived xenographs.
#'
#' @format A data frame with 38 rows and 3 variables:
#' \describe{
#'   \item{tumour}{The name of the PDX tumour}
#'   \item{ODC1}{The transcript levels of the ODC1 gene in the tumour}
#'   \item{angle}{The angle between the treatment and control}
#' }
#' @source Mer A, Haibe-Kains B (2021). Xeva: Analysis of patient-derived
#' xenograft (PDX) data. R package version 1.10.0.
"BRCA.PDXE.paxlitaxel.response"


#' Drug sensitivity data for genes in the GRCh38.p13 genome assembly
#'
#' A dataset containing drug sensitivity signatures for cancer cell-line
#' experiments. Each row contains the computed drug sensitivity signature
#' (estimate, p-value, fdr, n) of a single gene in a single experiment (an
#' experiment is defined by a tissue, compound/drug, and molecular data type).
#'
#' @format A data frame with 7805 rows and 12 variables:
#' \describe{
#'   \item{gene}{ENSEMBL gene name}
#'   \item{compound}{name of compound used in experiment}
#'   \item{tissue}{name of tissue used in experiment}
#'   \item{estimate}{The drug sensitivity estimate}
#'   \item{n}{The number of trials}
#'   \item{fdr}{false discovery rate (adjusted p-value)}
#'   \item{mDataType}{the molecular data type, i.e., rna, cnv, etc.}
#'   \item{gene_seq_start}{the relative starting position of the gene within
#'   the chromosome}
#'   \item{gene_seq_end}{the relative ending position of the gene within
#'   the chromosome}
#'   \item{strand}{which strand of DNA the gene is located on (+/-)}
#'   \item{chr}{the chromosome on which the gene is located}
#'   \item{pvalue}{the p-value of the drug response}
#' }
#' @source Feizi, N., Nair, S. K., Smirnov, P., Beri, G., Eeles, C., Esfahani,
#' P. N., … Haibe-Kains, B. (2021). PharmacoDB 2.0: Improving scalability and
#' transparency of in vitro pharmacogenomics analysis. _bioRxiv._
#' doi:10.1101/2021.09.21.461211
"Biomarkers"


#' GRCh38.p13 genome assembly
#'
#' A dataset containing the characteristics of of the GRCh38.p13 genome
#' assembly, specifically chromosome names and lengths.
#'
#' @format A data frame with 25 rows and 3 variables:
#' \describe{
#'   \item{chrName}{the name of the chromosome (all, 1-22, X, Y)}
#'   \item{molecule-type/loc}{the molecule type ("chromosome" or "all")}
#'   \item{chrLength}{the length of the chromosome (# nucleotides)}
#' }
#' @source Frankish, A., Diekhans, M., Ferreira, A. M., Johnson, R., Jungreis,
#' I. Loveland, J., Mudge, J. M., Sisu, C., Wright, J., Armstrong, J., Barnes,
#' I., Berry, A., Bignell, A., Carbonell Sala, S., Chrast, J., Cunningham, F.,
#' Di Domenico, T., Donaldson, S., Fiddes, I. T., García Girón, C., … Flicek, P.
#' (2019). GENCODE reference annotation for the human and mouse genomes.
#' _Nucleic acids research, 47_(D1), D766–D773.
#' https://doi.org/10.1093/nar/gky955
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
#' @source Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B.
#' L., Gillette, M. A., … Mesirov, J. P. (2005). Gene set enrichment analysis:
#' A knowledge-based approach for interpreting genome-wide expression profiles.
#' _Proceedings of the National Academy of Sciences, 102_(43), 15545–15550.
#' doi:10.1073/pnas.0506580102
"TestGeneSets"

# [END]
