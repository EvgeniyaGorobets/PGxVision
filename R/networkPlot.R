#' Build a network plot of gene sets
#'
#' Create a network plot where each node corresponds to a gene set and edge
#' weights correspond to the similarity score between two gene sets.
#'
#' @param gsSimilarityDf A data.frame listing the similarity score for each
#' pair of gene sets. Must have columns "gs1", "gs2", "similarity"
#' @param similarityCutoff (optional) A number indicating the minimum
#' similarity two gene sets must have in order for an edge to show up on the
#' plot. Defaults to 0.5. The higher this value is, the more legible the plot
#' will be.
#' @param title (optional) A custom title for the network plot. Defaults to
#' "Gene Set Similarity Plot"
#' @return invisible NULL
#'
#' @examples
#' geneSetIds <- queryGene("ENSG00000000971", "GO:BP")
#' geneSets <- expandGeneSets(geneSetIds, "GO:BP")
#' gsSimilarity <- computeGeneSetSimilarity(geneSets)
#' buildNetworkPlot(gsSimilarityDf, similarityCutoff=0.3)
#'
#' @importFrom checkmate assertDataFrame assertNames assertString
#' @importFrom igraph graph_from_data_frame E
#' @importFrom viridis magma
#' @export
buildNetworkPlot <- function(gsSimilarityDf, similarityCutoff=0.5, title=NULL) {
  # Check user input
  checkmate::assertDataFrame(gsSimilarityDf, min.rows=1)
  checkmate::assertNames(colnames(gsSimilarityDf),
                         must.include=c("gs1", "gs2", "similarity"))
  checkmate::assertNumber(similarityCutoff)

  if (!is.null(title)) {
    checkmate::assertString(title)
  } else {
    title <- "Gene Set Similarity Plot"
  }

  # Prune edges with similarity below the cutoff
  networkDf <- gsSimilarityDf[gsSimilarityDf$similarity >= similarityCutoff, ]
  if (nrow(networkDf) == 0) {
    stop(paste0("No edges were left in the network plot after applying a ",
                "similarity cutoff of ", similarityCutoff, ".\n",
                "Try a lower cutoff value (note that the max similarity in ",
                "this data.frame is ", max(gsSimilarityDf$similarity), ")."))
  }

  # Create network plot
  networkGraph <- igraph::graph_from_data_frame(networkDf, directed = F)
  igraph::E(networkGraph)$weight <- networkDf$similarity

  # Create color scale based on edge weights
  # Using magma color scale but removing the lighter (yellow) colors
  # for better edge visibility
  scaleColors <- colorRamp(rev(viridis::magma(10))[3:10])
  edgeRGBColors <- scaleColors(igraph::E(networkGraph)$weight)
  edgeHexColors <- apply(edgeRGBColors, MARGIN=1,
    FUN=function(x) {rgb(x[1], x[2], x[3], maxColorValue=255)})
  igraph::E(networkGraph)$color <- edgeHexColors

  # Convert to a plot
  plot(networkGraph, arrow.mode="-", vertex.label.color="black", vertex.size=20,
       edge.width=igraph::E(networkGraph)$weight*10,
       edge.color=igraph::E(networkGraph)$color)
  title(title)
  return(invisible(NULL)) #TODO: return something more useful
}


#' Perform a gene set analysis on a gene
#'
#' Query MSigDb to find all gene sets that queryGene is in, then compute the
#' similarity of the gene sets and plot them on a network graph. This is a
#' wrapper function which performs the full gene set analysis pipeline and
#' plots the result.
#'
#' @param geneId The ENSEMBL ID of the gene you want to query
#' @param queryType The type of query you want to perform; see MSigDb for
#' possible subcategories (use msigdbr::msigdbr_collections())
#' @param similarityCutoff (optional) A number indicating the minimum
#' similarity two gene sets must have in order for an edge to show up on the
#' plot. Defaults to 0.5.
#' @param title (optional) A custom title for the network plot. Defaults to
#' "querytype Gene Sets containing geneId (edge weights based on gene set
#' overlap)".
#' @return invisible(NULL)
#'
#' @examples
#' geneSetAnalysis("ENSG00000000971", "GO:BP", 0.3)
#'
#' @importFrom checkmate assertDataFrame assertNames assertString
#' @importFrom igraph graph_from_data_frame E
#' @importFrom viridis magma
#' @export
geneSetAnalysis <- function(geneId, queryType, similarityCutoff=NULL,
                            title=NULL) {
  # Perform gene set analysis
  geneSetIds <- queryGene(geneId, queryType)
  if (length(geneSetIds) <= 1) {
    stop(paste0("1 or fewer gene sets of type ", queryType, " and containing ",
                geneId, " were found.\n",
                "Further gene set analysis will not possible.\n",
                "Try a different gene or a different query type."))
  }

  geneSets <- expandGeneSets(geneSetIds, queryType)
  gsSimilarity <- computeGeneSetSimilarity(geneSets)

  if (is.null(title)) {
    title <- paste0(queryType, " Gene Sets containing ", geneId,
                    "\n(edge weights based on gene set overlap)")
  }

  buildNetworkPlot(gsSimilarity, similarityCutoff=similarityCutoff, title=title)
  return(invisible(NULL))
}
