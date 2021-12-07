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
#'
#' @return The visNetwork object containing the network plot of gene sets
#'
#' @examples
#' geneSetIds <- queryGene("ENSG00000000971", "GO:BP")
#' geneSets <- expandGeneSets(geneSetIds, "GO:BP")
#' gsSimilarityDf <- computeGeneSetSimilarity(geneSets)
#' buildNetworkPlot(gsSimilarityDf, similarityCutoff=0.3)
#'
#' @importFrom checkmate assertDataFrame assertNames assertString
#' @importFrom visNetwork visNetwork visPhysics visEvents
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

  # Create nodes
  geneSets <- unique(c(gsSimilarityDf$gs1, gsSimilarityDf$gs2))
  nodes <- data.frame(id = geneSets, label = geneSets, shape = "circle",
                      color="orange", font = "24px times black")

  # Remove edges between a node & itself and edges with similarity below cutoff
  networkDf <- gsSimilarityDf[gsSimilarityDf$gs1 != gsSimilarityDf$gs2,]
  networkDf <- networkDf[networkDf$similarity >= similarityCutoff, ]
  if (nrow(networkDf) == 0 && length(geneSets) > 1) {
    warning(paste0("No edges were left in the network plot after applying a ",
                   "similarity cutoff of ", similarityCutoff, ".\n",
                   "Try a lower cutoff value (note that the max similarity in ",
                   "this data.frame is ", max(gsSimilarityDf$similarity), ")."))
  }

  # Create color palette for edges
  scaleColors <- colorRamp(rev(viridis::magma(10))[3:10])
  edgeRGBColors <- scaleColors(networkDf$similarity)
  edgeHexColors <- apply(edgeRGBColors, MARGIN=1,
    FUN=function(x) {rgb(x[1], x[2], x[3], maxColorValue=255)})

  # Create edges
  edges <- data.frame(from = networkDf$gs1, to = networkDf$gs2,
                      value = networkDf$similarity, color = edgeHexColors)
  if (nrow(edges) > 0) {
    edges["smooth"] = F
  }

  # Build network
  network <- visNetwork::visNetwork(nodes, edges, main = title) %>%
    # Use physics to create proper layout
    visNetwork::visPhysics(
      solver = "forceAtlas2Based",
      forceAtlas2Based = list(gravitationalConstant = -255, avoidOverlap=1)) %>%
    # Once network has stabilized, disable physics so that users can drag
    # nodes around
    # Used JS code by Perry and YakovL:
    # Perry & YakovL. (2017). Stop vis.js physics after nodes load but allow
    # drag-able nodes. StackOverflow.
    # https://stackoverflow.com/questions/32403578/stop-vis-js-physics-after-nodes-load-but-allow-drag-able-nodes
    visNetwork::visEvents(stabilizationIterationsDone = "function () {
      network.setOptions( { physics: false } );}")

  return(network)
}

# [END]
