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
#' @return `ggplot` Plot object, showing the percentile of `value` relative to
#'   the distribution of `feature` in `reference`.
#'
#' @importFrom ggplot2 geom_line geom_area labs
#' @importFrom checkmate checkNumeric
#' @export
densityPlotBiomarkerPercentile <- function(feature, value, reference) {
    percentile <- computeBiomarkerPercentile(feature, value, reference)
    feature_vector <- as.numeric(reference[feature, ])
    z_score <- quantile(feature_vector, percentile)
    df <- data.frame(density(feature_vector)[c("x", "y")])
    cutoff <- df$x <= z_score
    if (percentile > 0.5) cutoff <- !cutoff
    plot <- ggplot(df, aes(x, y)) +
        geom_area(data=df[cutoff, ], fill="red") +
        geom_line() +
        labs(
            title=paste0(feature, " is ", round(percentile * 100, 1),
                "th Percentile Relative to Population"),
            x="Z Score",
            y="Frequency"
        )
    return(plot)
}