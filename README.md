
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

data.table is used to transform all data.frames into data.tables, for
faster manipulation and cleaner syntax. ggplot2 and ggprism are used to
create all plots in this package. checkmate is used to check user input.

## References

TODO

## Acknowledgements

This package was developed as part of an assessment for 2021 BCB410H:
Applied Bioinformatics, University of Toronto, Toronto, CANADA.
