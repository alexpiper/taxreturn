% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{alignment_entropy}
\alias{alignment_entropy}
\title{Alignment entropy}
\usage{
alignment_entropy(
  x,
  mask_gaps = 0.2,
  count_gaps = FALSE,
  method = "ML",
  unit = "log",
  return_extra = FALSE
)
}
\arguments{
\item{x}{A DNAbin or AAbin object}

\item{mask_gaps}{The threshold of gaps allowed before a position in the alignment is masked}

\item{count_gaps}{Whether gaps should be counted within entropy calculations. Default is FALSE.}

\item{method}{the method employed by \code{entropy::entropy} to estimate alignment entropy. Accepts:
"ML" : maximum likelihood
"MM" : bias-corrected maximum likelihood,
"Jeffreys" : Dirichlet with a=1/2
"Laplace" : Dirichlet with a=1
"SG" : Dirichlet with a=a=1/length(y)
"minimax" : Dirichlet with a=sqrt(sum(y))/length(y)
"CS" : ChaoShen
"NSB": Nemenman, Shafee and Biale (2002)
"shrink" : Shrinkage estimator
See the help page of \code{entropy::entropy} for more information}

\item{unit}{the unit in which entropy is measured. The default is "nats" (natural units). For computing entropy in "bits" set unit="log2".}

\item{return_extra}{Whether to return a dataframe including extra columns including individual base counts, gap proportions and number of bases at each position}
}
\description{
Alignment entropy
}
