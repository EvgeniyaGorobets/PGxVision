#' Build a volcano plot of biomarkers
#'
#' Note that this volcano plot is NOT plotting results of differential gene
#' expression analysis but is plotting the drug sensitivity signature (drug
#' sensitivity estimate on x-axis and log10(p-value) on y-axis). Points are
#' colored based on significance (red = significant; gray = not significant).
#' Currently, only a single experiment (one tissue, one compound, one molecular
#' data type) can be plotted at once.
#'
#' @param biomarkerDf A data.frame of drug sensitivity biomarkers that you want
#' to plot, including the compound, tissue, molecular data type ("mDataType"),
#' drug sensitivity estimate ("estimate") and p-value ("pvalue")
#' @param tissue A string representing the tissue for which you want to plot
#' biomarker drug sensitivity signatures
#' @param compound A string representing the compound for which you want to plot
#' biomarker drug sensitivity signatures
#' @param mDataType A string representing the molecular data type for which you
#' want to plot biomarker drug sensitivity signatures
#' @param pValCutoff (optional) A number between 0 and 1 indicating the maximum
#' p-value that an experimental result can have to be considered significant.
#' Any results with a higher p-value will be grayed out on the plot. Default
#' value is 0.05.
#' @param xLabel (optional) The label for the x-axis. Defaults to "Estimate".
#' @param yLabel (optional) The label for the y-axis. Defaults to
#' "Log10 P-Value".
#' @param title (optional) The title for the plot. Defaults to "Drug
#' sensitivity in <tissueName> tissue in response to <compoundName>
#'
#' @return A named list with two fields: "dt" and "plot". The "plot" field
#' contains a ggplot2 plot object (the volcano plot) and the "dt" field
#' contains a data.table with all the biomarkers that were plotted on "plot"
#'
#' @examples
#' data(Biomarkers)
#' result <- buildVolcanoPlot(Biomarkers, "Lung", "Trametinib", "rna", 0.005)
#' result$dt
#' result$plot
#'
#' @importFrom data.table setDT copy :=
#' @importFrom checkmate assertDataFrame assertNames assertNumber assertString
#' @importFrom ggplot2 ggplot geom_point scale_x_continuous theme aes
#' scale_color_manual ggtitle element_text geom_hline ylab xlab
#' @export
buildVolcanoPlot <- function(biomarkerDf, tissue="", compound="", mDataType="",
                             pValCutoff=0.05, xLabel=NULL, yLabel=NULL,
                             title=NULL) {
  # Local bindings to satisfy check() and DT syntax
  significant <- pvalue <- estimate <- NULL

  # Assign axis labels and title, if user did not provide any
  if (is.null(xLabel)) {
    xLabel <- "Estimate"
  }
  if (is.null(yLabel)) {
    yLabel <- "Log10 P-Value"
  }
  if (is.null(title)) {
    title <- paste("Drug sensitivity in", tissue, "tissue in response to",
                   compound)
  }

  # Check user inputs
  checkmate::assertDataFrame(biomarkerDf)
  checkmate::assertNames(colnames(biomarkerDf), must.include=c("tissue",
    "compound", "mDataType", "pvalue", "estimate"))
  checkmate::assertNumber(pValCutoff, lower=0, upper=1)
  checkmate::assertString(xLabel)
  checkmate::assertString(yLabel)
  checkmate::assertString(title)

  # Convert biomarkerDf to data.table and extract relevant biomarkers
  data.table::setDT(biomarkerDf, keep.rownames=TRUE)
  selectedBiomrkrs <- selectExperiment(biomarkerDf, tissue, compound, mDataType)

  # Create a copy and add a column which indicates whether each result is
  # significant, based on the p-value cutoff given
  selectedBiomrkrs <- data.table::copy(selectedBiomrkrs)
  selectedBiomrkrs[, significant := (pvalue <= pValCutoff)]
  selectedBiomrkrs[, log10pValue := -log10(pvalue)] #FIXME: these can be combined into one statement

  # Build the volcano plot
  plot <- ggplot2::ggplot(selectedBiomrkrs, ggplot2::aes(
    x=estimate, y=log10pValue, col=significant)) +
    # Add horizontal line to show significance cutoff
    ggplot2::geom_hline(yintercept=-log10(pValCutoff), linetype='dotted',
                        col = 'black', size=1) +
    # Add scatter points and color appropriately
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values=c("gray", "red")) +
    # Add title and axis labels and format them
    ggplot2::ggtitle(title) + ggplot2::ylab(yLabel) + ggplot2::xlab(xLabel) +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(hjust = 0.5))

  result <- list("dt" = selectedBiomrkrs, "plot" = plot)
  return(result)
}

# [END]
