#' Build a waterfall plot of ???
#'
#' More description???
#'
#' @param pdxDf A data.frame of ??
#' @return A ggplot2 plot object mapping ???
#'
#' @examples
#' data(PDXE)
#' df <- processXevaSet(PDXE, "BRCA", "paclitaxel", "RNASeq")
#' buildWaterfallPlot(df)
#'
#' @importFrom checkmate assertDataFrame assertNames
#' @importFrom ggplot2 ggplot geom_bar scale_fill_continuous theme aes
#' theme_classic ggtitle element_text ylab
#' @export
buildWaterfallPlot <- function(pdxDf) {
  # TODO: this already has selected experiment but maybe i should implement
  # so you can still select experiment

  # Check user inputs
  checkmate::assertDataFrame(pdxDf)
  checkmate::assertNames(colnames(pdxDf),
                         must.include=c("tumour", "ODC1", "angle"))

  # Build the waterfall (bar) plot
  plot <- ggplot(pdxDf, aes(x=tumour, y=angle, fill=ODC1))
  plot <- plot + geom_bar(stat="identity") +
          scale_fill_continuous(type="viridis")

  # Add title and axes labels
  # TODO: compound name shouldn't be hardcoded; more descriptive title?
  plot <- plot + theme_classic() + ggtitle("Paclitaxel Response") +
          ylab("Angle Between Treatment and Control") +
          theme(axis.text.x=element_blank(),
                plot.title = element_text(hjust = 0.5))

  return(plot)
}



# THis whole function is a mess and I should probably delete it
#' Prepare a XevaSet for input into buildWaterfallPlot
#'
#' @examples
#' data(PDXE)
#' processXevaSet(PDXE, "BRCA", "paclitaxel", "RNASeq")
#'
#' @importFrom Xeva subsetXeva summarizeData
#' @importFrom Biobase exprs
processXevaSet <- function(xevaSet, tissue, compound, mDataType) {
  # TODO: check xeva set somehow?
  # TODO: RENAME VARIABLES

  # Subset XevaSet
  brca <- Xeva::subsetXeva(PDXE, tissue, "tissue")
  # Retrieve molecular profiles
  brca_exprs <- Xeva::summarizeData(brca, drug=compound, mDataType=mDataType)

  # TODO: what is bcra_exprs$angle?
  # Order colnames by angle??
  orderedCols <- colnames(brca_exprs)[order(brca_exprs$angle)]
  # What is a facotor? Why is it called tumor? i don't understand anything
  tumour <- factor(orderedCols, levels=orderedCols)
  odc1 <- Biobase::exprs(brca_exprs)["ODC1",order(brca_exprs$angle)]
  x <- -brca_exprs$angle[order(brca_exprs$angle)]

  # Construct data.frame
  xevaDf <- data.frame(tumour=tumour, ODC1=odc1, angle=x, check.names=F)

  return(xevaDf)
}
