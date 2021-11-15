#' Build a network plot of gene sets
#'
#' Create a network plot where each node corresponds to a gene set and edge
#' weights correspond to the similarity score between two gene sets.
#'
#' @param gsSimilarityDf A data.frame listing the similarity score for each
#' pair of gene sets
#' @param title (optional) A custom title for the network plot. Defaults to
#' "Gene Set Similarity Plot"
#' @return
#'
#' @examples
#' geneSetNames <- queryGene("ENSG00000000971", "GO:BP")
#' geneSets <- expandGeneSets(geneSetNames, "GO:BP")
#' computeGeneSetSimilarity(geneSets)
#'
#' @importFrom checkmate assertDataFrame assertNames assertString
#' @importFrom igraph graph_from_data_frame E
#' @importFrom data.table setDT is.data.table
#' @export
buildNetworkPlot <- function(gsSimilarityDf, title=NULL) {
  # Chek user input
  checkmate::assertDataFrame(gsSimilarityDf, min.rows=1)
  checkmate::assertNames(colnames(gsSimilarityDf),
                         must.include=c("gs1", "gs2", "similarity"))
  if (!is.null(title)) {
    checkmate::assertString(title)
  } else {
    title <- "Gene Set Similarity Plot"
  }

  # Create network plot
  networkGraph <- igraph::graph_from_data_frame(gsSimilarityDf, directed = F)
  igraph::E(networkGraph)$weight <- gsSimilarityDf$similarity

  # Convert to a plot and return
  plot <- plot(networkGraph, arrow.mode="-",
               edge.width=igraph::E(networkGraph)$weight*10)
  return(plot)
}


temp <- function(geneId, queryType, title=NULL) {
  geneSetNames <- queryGene(geneId, queryType)
  geneSets <- expandGeneSets(geneSetNames, queryType)
  gsSimilarity <- computeGeneSetSimilarity(geneSets)

  if (is.null(title)) {
    title <- paste0(queryType, " Gene Sets containing ", geneId,
                    "\n(edge weights based on gene set overlap)")
  }
  plot <- buildNetworkPlot(gsSimilarity, title=title)
}
