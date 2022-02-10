#' Compute the percentile of a biomarker relative to some reference population.
#'
#' @param feature `character(1)` Feature name of the biomarkers of interest. Must
#'   be a valid row.name in the reference `data.frame`.
#' @param value `numeric(1)` Value of the feature, must be in the same units
#'   as the data in `reference`.
#' @param reference `data.frame` or `matrix` where rows are features, columns
#'   are samples, and values are a measurement on each feature for each sample.
#'   Rownames must match the specified `feature` value.
#'
#' @return `numeric(1)` Percentile of `value` relative to the distrubtion of
#'   `feature` in the `reference` population. Computed using `stats::ecdf`.
#'
#' @seealso `stats::ecdf`
#'
#' @importFrom checkmate checkCharacter checkNumeric checkDataFrame
#' @export
computeBiomarkerPercentile <- function(feature, value, reference) {
    suppressWarnings({
        reference_ecdf <- ecdf(reference[feature, ])
    })
    return(reference_ecdf(value))
}