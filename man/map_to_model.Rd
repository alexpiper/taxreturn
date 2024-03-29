% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_sequences.R
\name{map_to_model}
\alias{map_to_model}
\title{Map to model}
\usage{
map_to_model(
  x,
  model,
  min_score = 100,
  min_length = 1,
  max_indel = 9,
  max_gap = Inf,
  max_N = Inf,
  check_frame = FALSE,
  kmer_threshold = 0.5,
  k = 5,
  shave = TRUE,
  trim_ends = FALSE,
  extra = NA,
  multithread = FALSE,
  quiet = FALSE,
  progress = FALSE
)
}
\arguments{
\item{x}{A DNAbin or DNAStringset object}

\item{model}{A Profile Hidden Markov model ("PHMM" object) generated by \code{aphid::derivePHMM} to align the sequences to.
An already derived model of COI can be loaded using \code{data("model", package="taxreturn")}}

\item{min_score}{The minimum specificity (log-odds score for the optimal alignment) between the query sequence and the PHMM model for the sequence to be retained.
see \code{?aphid::Viterbi} for more information about the alignment process.}

\item{min_length}{The minimum length of the match between the query and PHMM for a sequence to retained. Takes into account sequential matches, as well as any internal insertions or deletions below max_indel.}

\item{max_indel}{The maximum number of internal insertions or deletions within the sequence to allow.}

\item{max_gap}{The maximum number of gaps within the sequence to allow.}

\item{max_N}{The max number of ambiguous N bases allowed before a sequence is removed.}

\item{check_frame}{Whether sequences with insertions or deletions which arent in multiples of 3 should be removed from output. Useful for coding loci such as COI but should not be used for non-coding loci. Default is FALSE.}

\item{kmer_threshold}{the maximum kmer distance allowed from the reference model. If a sequence is further than this, it will be skipped from the slower Viterbi alignment. Default is 50\% (0.5)}

\item{k}{integer giving the k-mer size used to generate the input matrix for k-means clustering. Default is k=5.}

\item{shave}{Whether bases that are outside (to the left or right) of the PHMM object should be removed from sequences in the output. Default is TRUE.}

\item{trim_ends}{Sometimes a trailing base can end up at the end of the alignment, separated by gaps. the trim_ends parameter checks up to n bases from each end of the alignment and if gaps are detected, any trailing bases will be removed.}

\item{extra}{How to handle insertions which were not part of the PHMM model. 'drop' will truncate all sequences to the shortest alignment length, while 'fill' will use gaps to pad all sequences out to the longest alignment length.}

\item{multithread}{Whether multithreading should be used, if TRUE the number of cores will be automatically detected (Maximum available cores - 1), or provide a numeric vector to manually set the number of cores to use. Default is FALSE (single thread)}

\item{quiet}{Whether progress should be printed to the console. Note that this will add additional runtime.}

\item{progress}{Whether a progress bar is displayed.}
}
\description{
This function alignes sequences to a Profile Hidden Markov Model (PHMM) using the Viterbi algorithm in order to retain only the target loci.
This function can also be used to extract smaller subregions out of longer sequences, for instance extracting the COI barcode region from mitochondrial genomes.
In order to reduce the number of sequences for alignment using the computationally expensive Viterbi algorithm, a rapid kmer distance screen is first conducted to remove any sequence too far diverged from the reference PHMM.
Similarly, all sequences that are more than twice the length of the reference PHMM model are broken into chunks of the same length as the model, and a rapid kmer screen conducted. The sequence is then subset to the most similar chunk and its two adjacent chunks prior to Viterbi alignment.
}
