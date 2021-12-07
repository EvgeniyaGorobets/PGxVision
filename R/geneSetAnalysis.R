#' Perform a gene set analysis on a gene
#'
#' Query MSigDb to find all gene sets that query gene (geneId) is in, then
#' compute the similarity of the gene sets and return the result in a data
#' frame. This is a wrapper function which performs the full gene set
#' analysis pipeline, minus the network plot.
#'
#' @param geneId The ENSEMBL ID of the gene you want to query
#' @param queryType The type of query you want to perform; see MSigDb for
#' possible subcategories (use msigdbr::msigdbr_collections())
#' @param similarityMetric The algorithm used to compute the similarity between
#' gene sets. The default (and only option currently) is "overlap".
#'
#' @return TODO
#'
#' @examples
#' geneSetAnalysis("ENSG00000000971", "GO:BP")
#'
#' @importFrom checkmate assertDataFrame assertNames assertString
#' @export
geneSetAnalysis <- function(geneId, queryType, similarityMetric="overlap") {
  # TODO: user input check
  # Perform gene set analysis
  targetGeneSets <- getGeneSets(geneId, queryType)
  if (nrow(targetGeneSets) == 0) {
    stop(paste0("There were no gene sets in the ", queryType,
               " category that contain ", geneId,
               ". Try a different gene or a differnt gene set category."))
  }

  geneSets <- expandGeneSets(targetGeneSets$gs_id, queryType)
  gsSimilarity <- computeGeneSetSimilarity(geneSets, similarityMetric)
  return(list(geneSets = targetGeneSets, similarityDf = gsSimilarity))
}


#' Get gene sets based on a query gene
#'
#' Get all gene sets from MSigDb that contain the query gene. Gene sets can be
#' based on biological process, cellular component, pathways, etc.
#'
#' @param geneId The ENSEMBL ID of the gene you want to query
#' @param queryType The type of query you want to perform; see MSigDb for
#' possible subcategories (use msigdbr::msigdbr_collections())
#'
#' @return A data.table containing all the gene sets in the category queryType
#' and containing the query gene. data.table contains gene set ID, name, source,
#' URL, and description.
#'
#' @examples
#' getGeneSets("ENSG00000000971", "GO:BP")
#'
#' @importFrom checkmate assertString
#' @importFrom msigdbr msigdbr msigdbr_collections
#' @importFrom data.table setDT
#' @export
getGeneSets <- function(geneId, queryType) {
  # Local bindings to satisfy check() and DT syntax
  ensembl_gene <- NULL

  # Check user inputs
  checkmate::assertString(geneId)
  checkmate::assertString(queryType)
  allQueryTypes <- msigdbr::msigdbr_collections()$gs_subcat
  checkmate::assertSubset(queryType, choices=allQueryTypes)

  # Get all gene sets in the queryType category
  allGeneSets <- msigdbr::msigdbr(species = "Homo sapiens",
                                  subcategory = queryType)
  data.table::setDT(allGeneSets)

  # Keep gene sets containing the query gene
  targetGeneSets <- allGeneSets[ensembl_gene == geneId, ]
  return(targetGeneSets[, list(gs_id, gs_name, gs_exact_source, gs_url,
                               gs_description)])
}


#' Expand gene sets
#'
#' Get all genes that are in the MSigDb gene sets corresponding to geneSetNames.
#'
#' @param geneSetIds A character vector containing gene set IDs from MSigDb
#' @param geneSetType (optional) The types of gene sets in geneSetNames; see
#' MSigDb for possible subcategories (use msigdbr::msigdbr_collections()). The
#' geneSetNames alone are enough but providing the geneSetType may speed up
#' computation because less data will be retrieved from MSigDb.
#'
#' @return A data.table mapping all gene sets to their component genes (2
#' columns: "gs_id", "ensembl_gene")
#'
#' @examples
#' geneSets <- getGeneSets("ENSG00000000971", "GO:BP")
#' expandGeneSets(geneSets$gs_id, "GO:BP")
#'
#' @importFrom checkmate assertCharacter assertString assertSubset
#' @importFrom msigdbr msigdbr
#' @importFrom data.table setDT data.table rbindlist
#' @export
expandGeneSets <- function(geneSetIds, geneSetType=NULL) {
  # Local bindings to satisfy check() and DT syntax
  gs_id <- ensembl_gene <- NULL

  # Check user inputs
  checkmate::assertCharacter(geneSetIds, min.len = 1)
  if (!is.null(geneSetType)) {
    checkmate::assertString(geneSetType)
    allQueryTypes <- msigdbr::msigdbr_collections()$gs_subcat
    checkmate::assertSubset(geneSetType, choices=allQueryTypes)
  }

  # Get all genes from humans
  allGenes <- msigdbr::msigdbr(species = "Homo sapiens",
                               subcategory = geneSetType)
  data.table::setDT(allGenes)

  # Map each gene set name to the set of corresponding genes
  geneList <- vector(mode = "list", length = length(geneSetIds))
  for (i in seq_along(geneSetIds)) {
    geneList[[i]] <- allGenes[gs_id == geneSetIds[i],
                              list(gs_id, ensembl_gene)]
  }

  # Combine into a single data.table and return
  expandedGeneSets <- data.table::rbindlist(geneList, use.names=TRUE)
  return(expandedGeneSets)
}


#' Compute the similarity between gene sets
#'
#' @param geneSets A data.table or data.frame with columns "gs_id" (gene set
#' ID) and "ensembl_gene" (ENSEMBL gene ID) which lists the genes in each
#' gene set
#' @param similarityMetric (optional) The type of similarity metric to compute.
#' Currently, the only option is "overlap", which calculates the proportion of
#' intersecting genes to total genes between each pair of gene sets.
#'
#' @return A data.frame that lists the similarity score between each pair of
#' gene sets. There will be three columns: "gs1" (gene set 1), "gs2" (gene set
#' 2), and "similarity".
#'
#' @examples
#' targetGeneSets <- getGeneSets("ENSG00000000971", "GO:BP")
#' geneSets <- expandGeneSets(targetGeneSets$gs_id, "GO:BP")
#' computeGeneSetSimilarity(geneSets)
#'
#' @importFrom checkmate assertDataFrame
#' @importFrom msigdbr msigdbr
#' @importFrom data.table setDT is.data.table
#' @export
computeGeneSetSimilarity <- function(geneSets, similarityMetric="overlap") {
  # Check user inputs
  checkmate::assertDataFrame(geneSets, min.rows=1)
  checkmate::assertNames(colnames(geneSets), must.include=c("gs_id",
                                                            "ensembl_gene"))
  checkmate::assertString(similarityMetric, pattern="overlap")

  # Convert to data.table if data.frame given
  if (!data.table::is.data.table(geneSets)) {
    data.table::setDT(geneSets)
  }

  # Compute similarity score of each pair of gene sets
  geneSetIds <- unique(geneSets$gs_id)
  similarityDf <- data.frame()
  for (i in 1:length(geneSetIds)) {
    for (j in i:length(geneSetIds)) {
      gs1 <- geneSetIds[i]
      gs2 <- geneSetIds[j]
      similarity <- overlapDistance(geneSets, gs1, gs2)
      similarityDf <- rbind(similarityDf, data.frame("gs1"=gs1, "gs2"=gs2,
                                                     "similarity"=similarity))
    }
  }

  return(similarityDf)
}


#' Helper function for computeGeneSetSimilarity
#'
#' Compute the overlap distance between two gene sets, defined as (number of
#' overlapping genes) / (number of unique genes in the union of gene sets).
#'
#' @param geneSets A data.table with columns "gs_name" (gene set name) and
#' "ensembl_gene" (ENSEMBL gene ID) which lists the genes in each gene set
#' @param gs1 The name of gene set 1
#' @param gs2 The name of gene set 2
#'
#' @return The overlap distance, as defined in the description, between the two
#' gene sets
overlapDistance <- function(geneSets, gs1, gs2) {
  # Local bindings to satisfy check() and DT syntax
  gs_id <- ensembl_gene <- NULL

  # Get genes in each gene set
  geneSet1 <- geneSets[gs_id == gs1, ensembl_gene]
  geneSet2 <- geneSets[gs_id == gs2, ensembl_gene]

  # Compute overlap
  numOverlapping <- length(intersect(geneSet1, geneSet2))
  totalGenes <- length(unique(c(geneSet1, geneSet2)))
  overlap <- numOverlapping / totalGenes

  return(overlap)
}

# [END]
