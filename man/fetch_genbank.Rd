% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fetch_sequences.R
\name{fetch_genbank}
\alias{fetch_genbank}
\title{all \verb{[filter]} 	Retrieves everthing
Specified \verb{[property]} 	Formal binomial and trinomial
at or below species level \verb{[property]}
family \verb{[rank]} 	Rank-based query
taxonomy genome \verb{[filter]} 	Taxa with a direct link to a genome sequence
2009/10/21:2020 \verb{[date]} 	Date-bounded query
mammalia \verb{[subtree]} 	All taxa within the Mammalia
extinct \verb{[property]} 	Extinct organisms
Terminal \verb{[property]} 	Terminal nodes in the tree
loprovencyclife \verb{[filter]} 	Entries with LinkOut links to the Encyclopedia of Life}
\usage{
fetch_genbank(
  x,
  database = "nuccore",
  marker = c("COI[GENE]", "CO1[GENE]", "COX1[GENE]"),
  output = "h",
  min_length = 1,
  max_length = 2000,
  subsample = FALSE,
  chunk_size = 100,
  db = NULL,
  multithread = FALSE,
  quiet = FALSE,
  progress = FALSE,
  retry_attempt = 3,
  retry_wait = 5
)
}
\arguments{
\item{x}{A taxon name or vector of taxa to download sequences for}

\item{database}{The database to download from. For NCBI GenBank this currently onlt accepts the arguments 'nuccore' or 'genbank' which is an alias for nuccore.}

\item{marker}{The barcode marker used as a search term for the database.
If this is set to "mitochondria" or "mitochondrion" it will download full mitochondrial genomes. If set to "genome" it will download entire genomes only.}

\item{output}{The output format for the taxonomy in fasta headers.
Options include "h" for full heirarchial taxonomy (Accession;Domain;Phylum;Class;Order;Family;Genus;Species),
"binom" for just genus species binomials (Accession;Genus_species),
"bold" for BOLD taxonomic ID only (Accession;BoldTaxID),
"gb" for genbank taxonomic ID (Accession|GBTaxID),
"gb-binom" which outputs genbank taxonomic ID's and Genus species binomials, translating BOLD taxonomic ID's to genbank in the process (Accession|GBTaxID;Genus_species)
or "standard" which outputs the default format for each database. For bold this is \verb{sampleid|species name|markercode|genbankid}}

\item{min_length}{The minimum length of sequences to download}

\item{max_length}{The maximum length of sequences to download}

\item{subsample}{(Numeric) return a random subsample of sequences from the search.}

\item{chunk_size}{Split up the query into chunks of this size to avoid overloading API servers. if left NULL, the default will be 300}

\item{db}{a database file generated using \code{taxreturn::get_ncbi_taxonomy()}. Generated automatically if NULL.}

\item{multithread}{Whether multithreading should be used, if TRUE the number of cores will be automatically detected, or provided a numeric vector to manually set the number of cores to use}

\item{quiet}{Whether progress should be printed to the console.}

\item{progress}{A logical, for whether or not to print a progress bar when multithread is true. Note, this will slow down processing.}

\item{retry_attempt}{The number of query attempts in case of query failure due to poor internet connection.}

\item{retry_wait}{How long to wait between query attempts.}
}
\description{
Fetch sequences from genbank
}
