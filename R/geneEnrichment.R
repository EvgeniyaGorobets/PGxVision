#' Query a gene
#'
#' Get all gene sets from MSigDb that contain the query gene. Gene sets can be
#' based on biological process, cellular component, pathways, etc.
#'
#' @param geneId The ENSEMBL ID of the gene you want to query
#' @param queryType The type of query you want to perform; see MSigDb for
#' possible subcategories (use msigdbr::msigdbr_collections())
#' @return A character vector containing the names of all gene sets this gene
#' is part of.
#'
#' @examples
#' queryGene("ENSG00000000971", "GO:BP")
#'
#' @importFrom checkmate assertString
#' @importFrom msigdbr msigdbr msigdbr_collections
#' @importFrom data.table setDT
#' @export
queryGene <- function(geneId, queryType) {
  # Check user inputs
  checkmate::assertString(geneId)
  checkmate::assertString(queryType)
  allQueryTypes <- msigdbr::msigdbr_collections()$gs_subcat
  checkmate::assertSubset(queryType, choices=allQueryTypes)

  allGeneSets <- msigdbr::msigdbr(species = "Homo sapiens",
                                  subcategory = queryType)
  data.table::setDT(allGeneSets)
  targetGeneSets <- allGeneSets[ensembl_gene == geneId]
  geneSetNames <- targetGeneSets$gs_name # TODO: maybe use gs_id instead?

  return(geneSetNames)
}


#' Expand gene sets
#'
#' Get all genes that are in the MSigDb gene sets corresponding to geneSetNames.
#'
#' @param geneSetNames A character vector containing gene set names from MSigDb
#' @param geneSetType (optional) The types of gene sets in geneSetNames; see
#' MSigDb for possible subcategories (use msigdbr::msigdbr_collections()). The
#' geneSetNames alone are enough but providing the geneSetType may speed up
#' computation because less data will be retrieved from MSigDb.
#' @return A data.table mapping all gene sets to their component genes (2
#' columns: "gs_name", "ensembl_gene")
#'
#' @examples
#' geneSetNames <- queryGene("ENSG00000000971", "GO:BP")
#' expandGeneSets(geneSetNames, "GO:BP")
#'
#' @importFrom checkmate assertCharacter assertString assertSubset
#' @importFrom msigdbr msigdbr
#' @importFrom data.table setDT data.table rbindlist
#' @export
expandGeneSets <- function(geneSetNames, geneSetType=NULL) {
  # Check user inputs
  checkmate::assertCharacter(geneSetNames, min.len = 1)
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
  geneList <- vector(mode = "list", length = length(geneSetNames))
  for (i in seq_along(geneSetNames)) {
    geneList[[i]] <- allGenes[gs_name == geneSetNames[i],
                              .(gs_name, ensembl_gene)]
  }

  # Combine into a single data.table and return
  expandedGeneSets <- data.table::rbindlist(geneList, use.names=TRUE)
  return(expandedGeneSets)
}
