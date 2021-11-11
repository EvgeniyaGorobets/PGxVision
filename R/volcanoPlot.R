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
#' @return A ggplot2 plot object mapping the biomarkers of the experiment
#' (x-axis = genome position; y-axis = p-value or fdr)
#'
#' @examples
#' data(Biomarkers)
#' experiment <- c("Lung", "Trametinib", "rna")
#' names(experiment) <- c("tissue", "compound", "mDataType")
#' buildVolcanoPlot(Biomarkers, experiment)
#'
#' @importFrom data.table setDT
#' @importFrom ggplot2 ggplot geom_point scale_x_continuous guides theme aes scale_color_manual
#' @importFrom ggprism guide_prism_minor
#' @export
buildVolcacoPlot <- function(biomarkerDf, experiment) {
  #

  # TODO: check inputs
  # Convert biomarkerDf to data.table and extract relevant biomarkers
  setDT(biomarkerDf, keep.rownames=TRUE)
  selectedBiomrkrs <- selectExperiment(biomarkerDf, experiment)

  # Build the volcano plot
  plot <- ggplot(selectedBiomrkrs, aes(x=estimate, y=-log10(pvalue)))
  plot <- plot + geom_point()

  return(plot)
}
