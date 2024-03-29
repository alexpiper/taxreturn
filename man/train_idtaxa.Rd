% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tax_utilities.R
\name{train_idtaxa}
\alias{train_idtaxa}
\title{Train IDTAXA}
\usage{
train_idtaxa(
  x,
  max_group_size = 10,
  max_iterations = 3,
  allow_group_removal = TRUE,
  orient = FALSE,
  get_lineage = FALSE,
  db = NULL,
  quiet = FALSE
)
}
\arguments{
\item{x}{A DNAbin}

\item{max_group_size}{The maximum size of any taxonomic group. This can be set to Inf (infinity) to allow for an unlimited number of sequences per group.}

\item{max_iterations}{The number of iterations to conduct training for in order to identify  any training sequences whose assigned classifications completely disagree with their predicted classification.}

\item{allow_group_removal}{Whether sequences that are the last remaining representatives of an entire group in the training data can be removed.
This can occur if the entire group appears to be misplaced in the taxonomic tree.}

\item{orient}{Training sequences must all be in the same orientation. Set this to TRUE to reorient the sequences if you are unsure.}

\item{get_lineage}{Get full taxonomic lineage using reformat_hierarchy if not already present.}

\item{db}{A database file generated by \code{get_ncbi_taxonomy} or \code{get_ott_lineage}. Required if get_lineage is TRUE.}

\item{quiet}{Whether progress should be printed to console.}
}
\description{
Train IDTAXA
}
