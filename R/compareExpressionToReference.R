#' Compute the percentile of a biomarker relative to some reference population.
#'
#' @param feature `character(1)` Feature name of the biomarkers of interest. Must
#'   be a valid row.name in the reference `data.frame`.
#' @param value `numeric(1)` Value of the feature, must be in the same units
#'   as the data in `reference`.
#' @param reference `data.frame` or `matrix` where rows are features, columns
#'   are samples, and values are a measurement on each feature for each sample.
#'
#' @examples
#'
#' @importFrom checkmate checkCharacter checkNumeric checkDataFrame
computeBiomarkerPercentile <- function(feature, value, reference) {
    reference_ecdf <- ecdf(reference[feature, ])
    return(reference_ecdf(value))
}

#'
#'
#'
#'
rectangularToLongFormat <- function(object) {

}