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
#' buildWaterfallPlot(df, xAxisCol="tumour", drugSensitivityCol="angle",
#'                    colorCol="ODC1", xLabel="Tumour",
#'                    yLabel="Angle Between Treatment and Control",
#'                    title="Paclitaxel Response in BRCA Tumours")
#'
#' @importFrom checkmate assertDataFrame assertNames assertString
#' @importFrom ggplot2 ggplot geom_bar scale_fill_continuous theme aes
#' theme_classic ggtitle element_text ylab xlab
#' @export
buildWaterfallPlot <- function(biomarkerDf, xAxisCol, drugSensitivityCol,
                               colorCol=NULL, xLabel=NULL, yLabel=NULL,
                               title=NULL) {
  # TODO: this already has selected experiment but maybe i should implement
  # so you can still select experiment

  # Check user inputs
  checkmate::assertDataFrame(biomarkerDf)
  checkmate::assertString(xAxisCol)
  checkmate::assertString(drugSensitivityCol)
  requiredCols <- c(xAxisCol, drugSensitivityCol)
  if (checkmate::testString(colorCol)) {
    requiredCols[2] <- colorCol
  } else {
    message(paste0("You have not selected a color scheme for your waterfall ",
                   "plot.\nBars will be colored based on drug sensitivity."))
  }
  # Check that the dataframe actually has the data cols specified by the user
  checkmate::assertNames(colnames(biomarkerDf), must.include=requiredCols)

  # Assign axes labels, if needed
  if (is.null(xLabel)) {
    xLabel <- xAxisCol
    message(paste0("You have not provided a custom label for the x-axis.\n",
                   "The label will default to xAxisCol."))
  }
  if (is.null(yLabel)) {
    yLabel <- drugSensitivityCol
    message(paste0("You have not provided a custom label for the y-axis.\n",
                   "The label will default to drugSensitivityCol."))
  }
  if (is.null(title)) {
    title <- paste0(yLabel, " vs. ", xLabel)
    message(paste0("You have not provided a custom title for your plot.\n",
                   "The title will default to 'yLabel vs. xLabel'."))
  }

  # Order the x-axis data points based on their drug sensitivity
  sortedDf <- biomarkerDf[order(biomarkerDf[,drugSensitivityCol],
                                decreasing=TRUE), ]
  sortedXAxis <- sortedDf[, xAxisCol]

  # Build the waterfall (bar) plot
  if (is.null(colorCol)) {
    fill <- NULL
  } else {
    fill <- sortedDf[, colorCol]
  }
  plot <- ggplot(sortedDf, aes(x=factor(sortedXAxis, levels=sortedXAxis),
                               y=sortedDf[, drugSensitivityCol], fill=fill))
  plot <- plot + geom_bar(stat="identity") +
          scale_fill_continuous(type="viridis")

  # Add title and axes labels
  plot <- plot + theme_classic() + ggtitle(title) + ylab(yLabel) +
    xlab(xLabel) + theme(plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  return(plot)
}



# THis whole function is a mess and I should probably delete it
#' Prepare a XevaSet for input into buildWaterfallPlot
#'
#' @examples
#' data(PDXE)
#' processXevaSet(PDXE, "BRCA", "paclitaxel", "RNASeq")
#'
#' @importFrom Xeva subsetXeva summarizeMolecularProfiles
#' @importFrom Biobase exprs
processXevaSet <- function(xevaSet, tissue, compound, mDataType) {
  # TODO: check xeva set somehow?
  # TODO: RENAME VARIABLES

  # Subset XevaSet
  brca <- Xeva::subsetXeva(PDXE, tissue, "tissue")
  # Retrieve molecular profiles
  brca_exprs <- Xeva::summarizeMolecularProfiles(brca, drug=compound, mDataType=mDataType)

  # What is a factor? Why is it called tumor? i don't understand anything
  tumour <- colnames(brca_exprs)
  odc1 <- Biobase::exprs(brca_exprs)["ODC1",]
  x <- brca_exprs$slope #TODO: wrong

  # Construct data.frame
  xevaDf <- data.frame(tumour=tumour, ODC1=odc1, angle=x, check.names=F)

  return(xevaDf)
}
