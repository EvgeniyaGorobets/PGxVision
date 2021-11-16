
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

Lang, M. (2017). “checkmate: Fast Argument Checks for Defensive R
Programming.” *The R Journal 9*(1), 437-445.
<https://journal.r-project.org/archive/2017/RJ-2017-028/index.html>.

Dowle, M., and Srinivasan, A. (2021). data.table: Extension of
`data.frame`. R package version 1.14.2.
<https://CRAN.R-project.org/package=data.table>.

Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis.
Springer-Verlag New York. <https://ggplot2.tidyverse.org>.

Dawson, C. (2021). ggprism: A ‘ggplot2’ Extension Inspired by ‘GraphPad
Prism’. R package version 1.0.3.
<https://CRAN.R-project.org/package=ggprism>

Csardi, G., Nepusz, T. (2006). The igraph software package for complex
network research, *InterJournal, Complex Systems 1695.*
<https://igraph.org>

Dolgalev, I. (2021). msigdbr: MSigDB Gene Sets for Multiple Organisms in
a Tidy Data Format. R package version 7.4.1.
<https://CRAN.R-project.org/package=msigdbr>

Garnier, S., Ross, N., Rudis, R., Camargo, A. P., Sciaini, M., and
Scherer, C. (2021). Rvision - Colorblind-Friendly Color Maps for R. R
package version 0.6.2.

## Acknowledgements

This package was developed as part of an assessment for 2021 BCB410H:
Applied Bioinformatics, University of Toronto, Toronto, CANADA.
