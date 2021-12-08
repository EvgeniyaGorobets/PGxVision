#' Build a Manhattan plot of biomarkers
#'
#' Note that this Manhattan plot is NOT plotting a GWAS study. It is plotting
#' the significance of drug responses of potential biomarkers (y-axis =
#' log10(pvalue)). Biomarkers are plotted in order of their genomic coordinates
#' (x-axis = gene sequence start) and colored by chromosome. Points that are
#' significant are fully opaque; points that are not significant are
#' translucent. Currently, only a single experiment (one tissue, one compound,
#' one molecular data type) can be plotted at once.
#'
#' @param biomarkerDf A data.frame of drug sensitivity biomarkers that you want
#' to plot, including the compound, tissue, molecular data type ("mDataType"),
#' genomic coordinates and chromosome of each gene ("gene_seq_start", "chr"),
#' and p-values ("pvalue")
#' @param chromosomeDf A data.frame of the reference genome, including the
#' names and lengths of all chromosomes ("chrName", "chrLength"). The row order
#' should match the chromosome order in the reference genome. If there is an
#' "all" row (with the total genome length), it should be first.
#' @param tissue A string representing the tissue for which you want to plot
#' biomarker drug sensitivity
#' @param compound A string representing the compound for which you want to plot
#' biomarker drug sensitivity
#' @param mDataType A string representing the molecular data type for which you
#' want to plot biomarker drug sensitivity
#' @param pValCutoff (optional) A number between 0 and 1 indicating the maximum
#' p-value that an experimental result can have to be considered significant.
#' Any results with a higher p-value will be translucent (opacity = 0.5) on the
#' plot. Default value is 0.05.
#' @param relativeGenomeCoords (optional) A boolean that indicates whether the
#' genomic coordinates given in biomarkerDf are relative to the chromosome
#' (TRUE) or relative to the entire genome (absolute, FALSE). Default is TRUE.
#' @param xLabel (optional) The label for the x-axis. Defaults to "Estimate".
#' @param yLabel (optional) The label for the y-axis. Defaults to
#' "Log10 P-Value".
#' @param title (optional) The title for the plot. Defaults to "Significance of
#' drug response in <tissueName> tissue in response to <compoundName>"
#'
#' @return A named list with two fields: "dt" and "plot". The "plot" field
#' contains a ggplot2 plot object (the Manhattan plot) and the "dt" field
#' contains a data.table with all the biomarkers that were plotted on "plot"
#'
#' @examples
#' data(Biomarkers)
#' data(GRCh38.p13.Assembly)
#' result <- buildManhattanPlot(
#'   Biomarkers, GRCh38.p13.Assembly, "Lung", "Trametinib", "rna",
#'   pValCutoff=0.01, relativeGenomeCoords=TRUE, title="Sample Plot")
#' result$dt
#' result$plot
#'
#' @importFrom data.table setDT copy :=
#' @importFrom checkmate assertDataFrame assertNames assertNumber assertString
#' assertSubset
#' @importFrom ggplot2 ggplot geom_point scale_x_continuous guides theme aes
#' scale_color_manual ggtitle element_text geom_hline xlim
#' @importFrom ggprism guide_prism_minor
#' @importFrom grDevices rainbow
#' @export
buildManhattanPlot <- function(biomarkerDf, chromosomeDf,
                               tissue="", compound="", mDataType="",
                               pValCutoff=0.05, relativeGenomeCoords=TRUE,
                               xLabel=NULL, yLabel=NULL, title=NULL) {
  # Local bindings to satisfy check() and DT syntax
  pvalue <- seq_start <- abs_gene_seq_start <- significant <- chrLength <-
    chrName <- chr <- NULL

  # Assign axis labels and title, if user doesn't provide any
  if (is.null(xLabel)) {
    xLabel <- "Chromosome"
  }
  if (is.null(yLabel)) {
    yLabel <- "Log10 P-Value"
  }
  if (is.null(title)) {
    title <- paste("Significance of drug response in", tissue,
                   "\ntissue in response to", compound)
  }

  # Check user inputs
  checkmate::assertDataFrame(biomarkerDf)
  checkmate::assertNames(colnames(biomarkerDf),
    must.include=c("tissue", "compound", "mDataType", "pvalue",
                   "gene_seq_start", "chr"))
  checkmate::assertDataFrame(chromosomeDf)
  checkmate::assertNames(colnames(chromosomeDf),
                         must.include=c("chrName", "chrLength"))
  checkmate::assertNumber(pValCutoff, lower=0, upper=1)
  checkmate::assertSubset(relativeGenomeCoords, choices=c(TRUE, FALSE))
  checkmate::assertString(xLabel)
  checkmate::assertString(yLabel)
  checkmate::assertString(title)

  # Convert dfs to data.table by reference
  data.table::setDT(biomarkerDf, keep.rownames=TRUE)
  data.table::setDT(chromosomeDf, keep.rownames=TRUE)

  # Select drug sensitivity results from the experiment
  selectedBiomrks <- selectExperiment(biomarkerDf, tissue, compound, mDataType)

  # If there is an "all" row in the chromosomeDf, remove it
  if (chromosomeDf[1, chrName] == "all") {
    chromosomeDf <- chromosomeDf[2:nrow(chromosomeDf)]
  }

  # If the provided genomic coordinates are relative to the chromosome,
  # shift them to be relative to the genome and compute chromosome coords
  if (relativeGenomeCoords) {
    result <- absolutizeGenomicCoords(selectedBiomrks, chromosomeDf)
    selectedBiomrks <- result[[1]]
    chromosomeDf <- result[[2]]
  }

  # Add a column which indicates whether each result is significant, based on
  # the p-value cutoff given
  # TODO: cutoff should be <= or < ?
  selectedBiomrks[, significant := (pvalue <= pValCutoff)]
  selectedBiomrks[, log10pValue := -log10(pvalue)] #FIXME: these can be combined into one statement

  # Build the Manhattan plot
  totalGenomeLen <- chromosomeDf[nrow(chromosomeDf), seq_start + chrLength]
  plot <- ggplot2::ggplot(selectedBiomrks, ggplot2::aes(
    x=abs_gene_seq_start, y=log10pValue, color=chr, alpha=significant)) +
    # Add horizontal line to show significance cutoff
    ggplot2::geom_hline(yintercept=-log10(pValCutoff), linetype='dotted',
                        col = 'black', size=1) +
    # Add scatter points
    ggplot2::geom_point()

  # Add title and colors
  plot <- plot + ggplot2::ggtitle(title) +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::scale_color_manual(values=rainbow(nrow(chromosomeDf)))

  # Customize x-axis ticks and labels
  midChromosome <- (chromosomeDf$seq_start + chromosomeDf$seq_end)/2
  chromosomeNames <- chromosomeDf$chrName
  plot <- plot + ggplot2::scale_x_continuous(
    "Chromosome", breaks=midChromosome, minor_breaks=c(1, chromosomeDf$seq_end),
    labels=chromosomeNames, limits=c(1, totalGenomeLen),
    guide = guide_prism_minor())

  result <- list("dt" = selectedBiomrks, "plot" = plot)
  return(result)
}

# [END]
