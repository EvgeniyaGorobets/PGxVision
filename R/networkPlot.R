#' Build a network plot of gene sets
#'
#' Create a network plot where each node corresponds to a gene set and edge
#' weights and lengths correspond to the similarity score between two gene sets.
#' The purpose is to visualize which gene sets share the most genes. Edges
#' with low similarity scores can be pruned to help naturally divide gene sets
#' into groups.
#'
#' @param gsSimilarityDf A non-empty data.frame listing the similarity score
#' for each pair of gene sets (essentially a table of edges). Must have columns
#' "gs1", "gs2", "similarity". It is ok to provide edges between a node and
#' itself, e.g., in the case of network plots with a single node. Such edges
#' will be removed.
#' @param similarityCutoff (optional) A number between 0 and 1 indicating the
#' minimum similarity two gene sets must have in order for an edge to show up
#' on the plot. Defaults to 0.5. Very low values will cause a highly connected,
#' possibly confusing plot. Very high values may remove or mask natural gene set
#' groupings.
#' @param title (optional) A custom title for the network plot. Defaults to
#' "Gene Set Similarity Plot"
#'
#' @return An interactive visNetwork object containing the network plot of
#' gene sets.
#'
#' @examples
#' result <- geneSetAnalysis("ENSG00000000971", "GO:BP")
#' buildNetworkPlot(result$similarityDf, similarityCutoff = 0.3)
#'
#' @importFrom checkmate assertDataFrame assertNames assertString
#' @importFrom visNetwork visNetwork visPhysics visEvents
#' @importFrom viridis magma
#' @importFrom grDevices colorRamp rgb
#' @export
buildNetworkPlot <- function(gsSimilarityDf, similarityCutoff=0.5, title=NULL) {
  # Check user input
  checkmate::assertDataFrame(gsSimilarityDf, min.rows=1)
  checkmate::assertNames(colnames(gsSimilarityDf),
                         must.include=c("gs1", "gs2", "similarity"))
  checkmate::assertNumber(similarityCutoff, lower=0, upper=1)

  if (!is.null(title)) {
    checkmate::assertString(title)
  } else {
    title <- "Gene Set Similarity Plot"
  }

  # Create nodes
  geneSets <- unique(c(gsSimilarityDf$gs1, gsSimilarityDf$gs2))
  nodes <- data.frame(
    id = geneSets, label = geneSets, shape = "circle",
    color=list(background = "orange", border = "#CC5500",
               highlight = list(background = "#F7D527", border = "#EBBA18")),
    font = "24px times black", borderWidthSelected = 3)

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
      network.setOptions( { physics: false } );}")  %>%
    # For now, avoid dragging nodes because it is weird after we disable physics
    visInteraction(dragNodes = FALSE) #visOptions(highlightNearest = TRUE)

  return(network)
}

# [END]
