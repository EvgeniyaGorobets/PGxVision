
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PGxVision

<!-- badges: start -->
<!-- badges: end -->

## Description

PGxVision (PharmacoGenomic Vision & Interpretation) helps identify and
visualize RNA-based cancer biomarkers for drug response. This package is
intended to be used in conjunction with the Roche-PharmacoGx pipeline.
PGxVision is intended to guide cancer treatment decisions in molecular
tumour boards.

``` r
R version 4.1.1 (2021-08-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19043)
```

## Installation

To download the package:

``` r
# install.packages("devtools")
devtools::install_github("EvgeniyaGorobets/PGxVision")
library("PGxVision")
```

To run the shinyApp:

``` r
Under construction
```

## Overview

``` r
ls("package:PGxVision")
data(package = "PGxVision")
```

TODO: add descriptions! \* Manhattan Plot: \* Volcano Plot: \* Waterfall
Plot: TODO: add vignettes and include image

## Contributions

The author of the package is Evgeniya Gorobets.

data.table is used to transform data.frames into data.tables in some
plotting and gene set analysis functions for faster manipulation and
cleaner syntax (*buildVolcanoPlot, buildManhattanPlot, queryGene,
expandGeneSets, computeGeneSetSimilarity*).

ggplot2 is used to create non-network plots (buildVolcanoPlot,
buildManhattanPlot, buildWaterfallPlot), and igraph is used to create
network plots (buildNetworkPlot). ggprism is used to enhance the axes on
the Manhattan plot (buildManhattanPlot). viridis is used to enhance the
colors on the network plot (buildNetworkPlot).

checkmate is used to succinctly check user input in all functions.

msigdbr is used to query the MSigDb in gene set analysis functions
(queryGene, expandGeneSets).

## References

Csardi, G., Nepusz, T. (2006). The igraph software package for complex
network research, *InterJournal, Complex Systems 1695.*
<https://igraph.org>

Dawson, C. (2021). ggprism: A ‘ggplot2’ Extension Inspired by ‘GraphPad
Prism’. R package version 1.0.3.
<https://CRAN.R-project.org/package=ggprism>

Dolgalev, I. (2021). msigdbr: MSigDB Gene Sets for Multiple Organisms in
a Tidy Data Format. R package version 7.4.1.
<https://CRAN.R-project.org/package=msigdbr>

Dowle, M., and Srinivasan, A. (2021). data.table: Extension of
`data.frame`. R package version 1.14.2.
<https://CRAN.R-project.org/package=data.table>.

Garnier, S., Ross, N., Rudis, R., Camargo, A. P., Sciaini, M., and
Scherer, C. (2021). Rvision - Colorblind-Friendly Color Maps for R. R
package version 0.6.2.

Lang, M. (2017). “checkmate: Fast Argument Checks for Defensive R
Programming.” *The R Journal 9*(1), 437-445.
<https://journal.r-project.org/archive/2017/RJ-2017-028/index.html>.

R Core Team. (2021). R: A language and environment for statistical
computing. R Foundation for Statistical Computing, Vienna, Austria.
<https://www.R-project.org/>.

Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis.
Springer-Verlag New York. <https://ggplot2.tidyverse.org>.

## Acknowledgements

This package was developed as part of an assessment for 2021 BCB410H:
Applied Bioinformatics, University of Toronto, Toronto, CANADA.
