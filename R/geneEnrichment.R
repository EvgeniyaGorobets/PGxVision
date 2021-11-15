#' Query a gene
#'
#' Get all gene sets from MSigDb that contain the query gene. Gene sets can be
#' based on biological process, cellular  component, or molecular function.
#'
#' @param geneId The ENSEMBL ID of the gene you want to query
#' @param queryType The type of query you want to perform; either "GO:BP" (GO
#' biological pathway), "GO:CC" (GO cellular component), ...
#' @return A character vector containing the names of all gene sets this gene
#' is part of.
#'
#' @examples
#' queryGene("ENSG00000000971", "GO:BP")
#'
#' @importFrom checkmate assertString
#' @importFrom msigdbr msigdbr
#' @importFrom data.table setDT
#' @export
queryGene <- function(geneId, queryType) {
  # Check user inputs
  checkmate::assertString(geneId)
  checkmate::assertString(queryType, pattern="(GO:BP)|(GO:CC)")

  allGeneSets <- msigdbr::msigdbr(species = "Homo sapiens",
                                  subcategory = queryType)
  data.table::setDT(allGeneSets)
  targetGeneSets <- allGeneSets[ensembl_gene == geneId]
  geneSetNames <- targetGeneSets$gs_name # TODO: maybe use gs_id instead?

  return(geneSetNames)
}
