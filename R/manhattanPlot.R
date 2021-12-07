#' Build a Manhattan plot of biomarkers
#'
#' @param biomarkerDf A data.frame of drug sensitivity biomarkers that you want
#' to plot, including the compound, tissue, and molecular data type used in the
#' experiment, the genomic coordinates of the genes which were tested on, and
#' the resulting p-value/fdr
#' @param chromosomeDf A data.frame of the reference genome, including the
#' lengths of all the chromosomes. The row order should match the chromosome
#' order in the reference genome. If there is an "all" row, it should be first.
#' @param tissue A string representing the tissue for which you want to plot
#' biomarkers
#' @param compound A string representing the compound for which you want to plot
#' biomarkers
#' @param mDataType A string representing the molecular data type for which you
#' want to plot biomarkers
#' @param pValCutoff (optional) A decimal number indicating the cutoff value for
#' significant observations; any results with a higher p-value will be grayed
#' out on the plot. Default value is 0.05.
#' @param relativeGenomeCoords (optional) A boolean that indicates whether the
#' genomic coordinates given in biomarkerDf are relative to the chromosome
#' (TRUE) or relative to the entire genome (absolute, FALSE). Default is TRUE.
#' @param xLabel (optional) The label for the x-axis. Defaults to "Estimate".
#' @param yLabel (optional) The label for the y-axis. Defaults to
#' "Log10 P-Value".
#' @param title (optional) The title for the plot. Defaults to "Drug
#' sensitivity in <tissueName> tissue in response to <compoundName>
#'
#' @return A ggplot2 plot object mapping the biomarkers of the experiment
#' (x-axis = genome position; y-axis = -log10(p-value or fdr))
#' TODO: update!
#'
#' @examples
#' data(Biomarkers)
#' data(GRCh38.p13.Assembly)
#' buildManhattanPlot(Biomarkers, GRCh38.p13.Assembly, "Lung", "Trametinib",
#'                    "rna", pValCutoff=0.01, relativeGenomeCoords=TRUE)
#'
#' @importFrom data.table setDT copy :=
#' @importFrom checkmate assertDataFrame assertNames assertNumber
#' @importFrom ggplot2 ggplot geom_point scale_x_continuous guides theme aes
#' scale_color_manual ggtitle element_text geom_hline xlim
#' @importFrom ggprism guide_prism_minor
#' @importFrom grDevices rainbow
#' @export
buildManhattanPlot <- function(biomarkerDf=NULL, chromosomeDf=NULL,
                               tissue="", compound="", mDataType="",
                               pValCutoff=0.05, relativeGenomeCoords=TRUE,
                               xLabel=NULL, yLabel=NULL, title=NULL) {
  # Local bindings to satisfy check() and DT syntax
  pvalue <- seq_start <- abs_gene_seq_start <- significant <- chrLength <-
    chrName <- chr <- NULL

  # Check user inputs
  checkmate::assertDataFrame(biomarkerDf)
  checkmate::assertNames(colnames(biomarkerDf),
    must.include=c("tissue", "compound", "mDataType", "pvalue",
                   "gene_seq_start", "chr"))
  checkmate::assertDataFrame(chromosomeDf)
  checkmate::assertNames(colnames(chromosomeDf),
                         must.include=c("chrName", "chrLength"))
  checkmate::assertNumber(pValCutoff, lower=0, upper=1)

  # Assign axis labels and title, if needed
  if (is.null(xLabel)) {
    xLabel <- "Chromosome"
  }
  if (is.null(yLabel)) {
    yLabel <- "Log10 P-Value"
  }
  if (is.null(title)) {
    title <- paste("Significance of drug response in", tissue,
                   "tissue in response to", compound)
  }

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
    labels=chromosomeNames, limits=c(1, totalGenomeLen), guide = guide_prism_minor()) #+ ggplot2::guides(x = guide_prism_minor())

  # Add horizontal line to show significance cutoff
  plot <- plot + ggplot2::geom_hline(yintercept=-log10(pValCutoff),
                                     linetype='dotted', col = 'black', size=1)

  result <- list("dt" = selectedBiomrks, "plot" = plot)
  return(result)
}


absolutizeGenomicCoords <- function(selectedBiomrks, chromosomeDf) {
  # Local bindings to satisfy check() and DT syntax
  seq_end <- seq_start <- abs_gene_seq_start <- gene_seq_start <- chrLength <-
    chrName <- chr <- NULL

  # Make a copy of the data.tables so you're not modifying user data
  selectedBiomrks <- data.table::copy(selectedBiomrks)
  chromosomeDf <- data.table::copy(chromosomeDf)

  # Add columns to track absolute position in genome (TODO: can maybe remove this)
  chromosomeDf$seq_start <- 1
  chromosomeDf$seq_end <- 1
  selectedBiomrks$abs_gene_seq_start <- 1

  # Shift genomic coords to be relative to genome instead of each chromosome
  chromEnd <- 0  # the end of the previous chromosome
  for (row in 1:nrow(chromosomeDf)) {
    # Set the absolute seq_start & seq_end of each chromosome
    chromosomeDf[row, seq_start := chromEnd + 1]
    chromosomeDf[row, seq_end := chromEnd + chrLength]

    # Update all biomarkers in that chromosome with absolute gene_seq_start
    chromName <- paste("chr", chromosomeDf[row, chrName], sep="") #TODO: not generic enough
    selectedBiomrks[chr == chromName,
                    abs_gene_seq_start := gene_seq_start + chromEnd]

    # Update chromEnd to end of current chromosome
    chromEnd <- chromosomeDf[row, seq_end]
  }

  # Return the modified dfs
  return(list(selectedBiomrks, chromosomeDf))
}

# [END]
