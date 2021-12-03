#' Build a volcano plot of biomarkers
#'
#' Note that this volcano plot is NOT plotting results of differential gene
#' expression analysis but is plotting the drug sensitivity signature
#'
#' @param biomarkerDf A data.frame of drug sensitivity biomarkers that you want
#' to plot, including the compound, tissue, and molecular data type used in the
#' experiment, the genomic coordinates of the genes which were tested on, and
#' the resulting p-value/fdr and estimate
#' @param experiment A named character vector representing the experiment for
#' which you want to plot biomarkers; an experiment is defined by a "tissue",
#' "compound", and "mDataType" (molecular data type)
#' @param pValCutoff A decimal number indicating the cutoff value for
#' significant observations; any results with a higher p-value will be grayed
#' out on the plot. Default value is 0.05.
#' @return A ggplot2 plot object mapping the biomarkers of the experiment
#' (x-axis = estimate; y-axis = -log10(p-value))
#' #TODO: update
#'
#' @examples
#' data(Biomarkers)
#' experiment <- c("Lung", "Trametinib", "rna")
#' names(experiment) <- c("tissue", "compound", "mDataType")
#' buildVolcanoPlot(Biomarkers, experiment, 0.005)
#'
#' @importFrom data.table setDT copy :=
#' @importFrom checkmate assertDataFrame assertNames assertNumber
#' @importFrom ggplot2 ggplot geom_point scale_x_continuous theme aes
#' scale_color_manual ggtitle element_text geom_hline
#' @export
buildVolcanoPlot <- function(biomarkerDf, experiment, pValCutoff=0.05) {
  # Local bindings to satisfy check() and DT syntax
  significant <- pvalue <- estimate <- NULL

  # Check user inputs
  checkmate::assertDataFrame(biomarkerDf)
  checkmate::assertNames(colnames(biomarkerDf), must.include=c("tissue",
    "compound", "mDataType", "pvalue", "estimate"))
  checkmate::assertNumber(pValCutoff, lower=0, upper=1)

  # Convert biomarkerDf to data.table and extract relevant biomarkers
  data.table::setDT(biomarkerDf, keep.rownames=TRUE)
  selectedBiomrkrs <- selectExperiment(biomarkerDf, experiment)

  # Create a copy and add a column which indicates whether each result is
  # significant, based on the p-value cutoff given
  # TODO: cutoff should be <= or < ?
  selectedBiomrkrs <- data.table::copy(selectedBiomrkrs)
  selectedBiomrkrs[, significant := (pvalue <= pValCutoff)]
  selectedBiomrkrs[, log10pValue := -log10(pvalue)] #FIXME: these can be combined into one statement

  # Build the volcano plot
  plot <- ggplot2::ggplot(selectedBiomrkrs, ggplot2::aes(
    x=estimate, y=log10pValue, col=significant))
  plot <- plot + ggplot2::geom_point() +
    ggplot2::scale_color_manual(values=c("gray", "red")) +
    ggplot2::ggtitle(paste0("P-value vs. estimate of biomarkers in ",
                            experiment["tissue"], " tissue in response to ",
                            experiment["compound"])) +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::geom_hline(yintercept=-log10(pValCutoff), linetype='dotted',
                        col = 'black', size=1)

  result <- list("dt" = selectedBiomrkrs, "plot" = plot)
  return(result)
}

# [END]
