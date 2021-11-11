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
#' out on the plot. Default value is 0.01.
#' @return A ggplot2 plot object mapping the biomarkers of the experiment
#' (x-axis = estimate; y-axis = -log10(p-value))
#'
#' @examples
#' data(Biomarkers)
#' experiment <- c("Lung", "Trametinib", "rna")
#' names(experiment) <- c("tissue", "compound", "mDataType")
#' buildVolcanoPlot(Biomarkers, experiment, 0.005)
#'
#' @importFrom data.table setDT
#' @importFrom ggplot2 ggplot geom_point scale_x_continuous guides theme aes scale_color_manual
#' @importFrom ggprism guide_prism_minor
#' @export
buildVolcanoPlot <- function(biomarkerDf, experiment, pValCutoff=0.01) {
  # TODO: check inputs

  # Convert biomarkerDf to data.table and extract relevant biomarkers
  setDT(biomarkerDf, keep.rownames=TRUE)
  selectedBiomrkrs <- selectExperiment(biomarkerDf, experiment)

  # Create a copy and add a column which indicates whether each result is
  # significant, based on the p-value cutoff given
  # TODO: cutoff should be <= or < ?
  selectedBiomrkrs <- copy(selectedBiomrkrs)
  selectedBiomrkrs[, significant := (pvalue <= pValCutoff)]

  # Build the volcano plot
  plot <- ggplot(selectedBiomrkrs, aes(x=estimate, y=-log10(pvalue),
                                       col=significant))
  plot <- plot + geom_point() + scale_color_manual(values=c("gray", "red"))

  return(plot)
}
