% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/networkPlot.R
\name{buildNetworkPlot}
\alias{buildNetworkPlot}
\title{Build a network plot of gene sets}
\usage{
buildNetworkPlot(gsSimilarityDf, similarityCutoff = 0.5, title = NULL)
}
\arguments{
\item{gsSimilarityDf}{A non-empty data.frame listing the similarity score
for each pair of gene sets (essentially a table of edges). Must have columns
"gs1", "gs2", "similarity". It is ok to provide edges between a node and
itself, e.g., in the case of network plots with a single node. Such edges
will be removed.}

\item{similarityCutoff}{(optional) A number between 0 and 1 indicating the
minimum similarity two gene sets must have in order for an edge to show up
on the plot. Defaults to 0.5. Very low values will cause a highly connected,
possibly confusing plot. Very high values may remove or mask natural gene set
groupings.}

\item{title}{(optional) A custom title for the network plot. Defaults to
"Gene Set Similarity Plot"}
}
\value{
An interactive visNetwork object containing the network plot of
gene sets.
}
\description{
Create a network plot where each node corresponds to a gene set and edge
weights and lengths correspond to the similarity score between two gene sets.
The purpose is to visualize which gene sets share the most genes. Edges
with low similarity scores can be pruned to help naturally divide gene sets
into groups.
}
\examples{
result <- geneSetAnalysis("ENSG00000000971", "GO:BP")
buildNetworkPlot(result$similarityDf, similarityCutoff = 0.3)

}
