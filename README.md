
<!-- README.md is generated from README.Rmd. Please edit that file -->
taxreturn
=========

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/alexpiper/taxreturn.svg?branch=master)](https://travis-ci.org/alexpiper/taxreturn) <!-- badges: end -->

taxreturn is an R package for fetching DNA barcode sequences and associated taxonomic annotations from public databases such as the Barcode of Life Database (BOLD) and NCBI GenBank, curating these sequences and formatting them into training sets compatible with popular taxonomic classifiers used for metabarcoding and marker gene analysis.

Installation
------------

This package is still in development and not yet available on CRAN. You can install development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
#devtools::install_github("alexpiper/taxreturn")
```

Current progress:
=================

Database creation
-----------------

-   fetchSeqs
-   cleanSeqs - use a profile hidden markov model to clean non-homologous sequences
-   format\_datbase(options="blast","RDP","q2","IDTAXA" etc))
-   add\_sequences(add local db)
-   Trim\_to\_primer\_regions - Using insect r package HMM's

Sample Sheet & Index switching
------------------------------

-   Samplesheet\_new(options=Miseq, hiseq etc)
-   Samplesheet\_allcombinations(expand out all combinations)
-   switchrate\_allcomb(calculate switchrate from expected vs observed)
-   switchrate\_control(calculate switchrate from control - allow selecting control)

Primer evaluation
-----------------

-   evaluate\_primerbind(use hmm similar to insect r package to locate bind site, then primer miner evaluation)
-   plot\_primer - Plot the above output\*

Plotting
--------

-   Mostly hanadled by phyloseq
-   summarise\_Taxa(samples=c()) - similar to current summarise taxa but allow subsettign to desired samples
-   Plot\_expectedobserved(allow comparison of fit)

Summarising
-----------

-   output detection tables etc - reports
-   How many sequences were removed at each stage - Plot of this (Histogram of reductions)
-   How many sequences were dereplicated - Plot of dataset redundancy (Histogram of species had 1,2,3,4,5,6,7,8,9,10+ etc)

Examples
--------

This is a basic example which shows you how to solve a common problem
