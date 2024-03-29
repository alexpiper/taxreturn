% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{acc2hex}
\alias{acc2hex}
\title{Convert accession number to hexadecimal coding}
\usage{
acc2hex(x, force = FALSE)
}
\arguments{
\item{x}{A DNAbin or DNAStringSet with names formatted Accession|taxid;taxonomy.
Or a character vector with names formatted Accession|taxid;taxonomy
Or a character vector of accessions}

\item{force}{override checks if string is already hexadecimal}
}
\description{
Hexadecimal coding is useful for ensuring that accession numbers dont contain characters that are considered illegal for newick strings
}
