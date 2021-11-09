# For development only, load sample data
if (sys.nframe() == 0) {
  # development code here
  source("R/utils.R")

  # TODO: example data should be saved as an R obj? why?
  biomarkerFile <- "data/pharmacodb_biomarkers_by_tissue.csv"
  biomarkerDf <- read.csv(biomarkerFile)

  experiment <- c("Lung", "Trametinib", "rna")
  names(experiment) <- c("tissue", "compound", "mDataType")
}

#' @importFrom data.table setDT
#' @importFrom ggplot2 ggplot, geom_point
buildVolcacoPlot <- function(biomarkerDf, experiment) {
  # NOTE: not actually plotting results of differential analysis
  # Used to look at drug sensitivity signature computed by PharmacoGx

  # TODO: check inputs
  # Convert biomarkerDf to data.table and extract relevant biomarkers
  setDT(biomarkerDf, keep.rownames=TRUE)
  selectedBiomrkrs <- selectExperiment(biomarkerDf, experiment)

  # Build the volcano plot
  plot <- ggplot(selectedBiomrkrs, aes(x=estimate, y=-log10(pvalue)))
  plot <- plot + geom_point()

  return(plot)
}
