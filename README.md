
<!-- README.md is generated from README.Rmd. Please edit that file -->

# taxreturn

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/alexpiper/taxreturn.svg?branch=master)](https://travis-ci.org/alexpiper/taxreturn)
<!-- badges: end -->

taxreturn is an R package for fetching DNA barcode sequences and
associated taxonomic annotations from public databases such as the
Barcode of Life Database (BOLD) and NCBI GenBank, curating these
sequences and formatting them into training sets compatible with popular
taxonomic classifiers used for metabarcoding and marker gene analysis.

**This package is still in development and many functions are likely to
change**

## Installation

This package is still in development and not yet available on CRAN. You
can install development version from [GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("alexpiper/taxreturn")
library(taxreturn)
```

## Examples

The main vignette demonstrates example usage of all main functions. This
can be accessed using the following code, or a rendered version can be
accessed
[here](https://alexpiper.github.io/taxreturn/doc/taxreturn-vignette.html)

``` r
# Use this to view the vignette as an isolated HTML file
utils::browseVignettes(package = "taxreturn")
```
