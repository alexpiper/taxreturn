# QC Max EE plot  ----------------------------------------------------------

#' QC Max EE plot
#' @description modified from code by Remi Maglione https://github.com/RemiMaglione/r-scripts/blob/master/Qual_vs_MaxEE_plot.R
#' @param fastq.r1
#' @param fastq.r2
#' @param name1
#' @param name2
#'
#' @return
#'
#' @examples
qualMaxEEplot <- function(fastq.r1, fastq.r2, name1 = "Fastq_R1", name2 = "Fastq_R2") {
  ### Loading dependenies
  require(ggplot2)
  require(patchwork)
  ### Definning function
  fastqTmp <- function(fastq) {
    fastq.tmp <- rbind(
      data.frame(
        R = fastq$X.Base,
        Q = fastq$Mean, S = c("Mean"),
        E = 10^(-fastq$Mean / 10),
        A = Reduce("+", 10^(-fastq$Mean / 10), accumulate = TRUE)
      ),
      data.frame(
        R = fastq$X.Base,
        Q = fastq$Median,
        S = c("Median"),
        E = 10^(-fastq$Median / 10),
        A = Reduce("+", 10^(-fastq$Median / 10), accumulate = TRUE)
      ),
      data.frame(
        R = fastq$X.Base,
        Q = fastq$Lower.Quartile,
        S = c("Lower.Quartile"),
        E = 10^(-fastq$Lower.Quartile / 10),
        A = Reduce("+", 10^(-fastq$Lower.Quartile / 10), accumulate = TRUE)
      ),
      data.frame(
        R = fastq$X.Base,
        Q = fastq$Upper.Quartile,
        S = c("Upper.Quartile"),
        E = 10^(-fastq$Upper.Quartile / 10),
        A = Reduce("+", 10^(-fastq$Upper.Quartile / 10), accumulate = TRUE)
      ),
      data.frame(
        R = fastq$X.Base,
        Q = fastq$X10th.Percentile,
        S = c("X10th.Percentile"),
        E = 10^(-fastq$X10th.Percentile / 10),
        A = Reduce("+", 10^(-fastq$X10th.Percentile / 10), accumulate = TRUE)
      ),
      data.frame(
        R = fastq$X.Base,
        Q = fastq$X90th.Percentile,
        S = c("X90th.Percentile"),
        E = 10^(-fastq$X90th.Percentile / 10),
        A = Reduce("+", 10^(-fastq$X90th.Percentile / 10), accumulate = TRUE)
      )
    )
    return(fastq.tmp)
  }

  qualPlot <- function(df.tmp) {
    p_r <- ggplot(df.tmp, aes(color = S)) +
      geom_point(aes(x = R, y = Q), size = 1) +
      labs(x = "Reads position", y = "Reads Quality")
    return(p_r)
  }

  maxEEplot <- function(df.tmp) {
    q_r <- ggplot(df.tmp[complete.cases(df.tmp), ], aes(color = S)) +
      geom_point(aes(x = R, y = log10(A)), size = 1) +
      geom_hline(yintercept = log10(2), color = "red") +
      geom_hline(yintercept = log10(3), color = "red") +
      geom_hline(yintercept = log10(5), color = "red") +
      geom_hline(yintercept = log10(7), color = "red") +
      geom_text(label = "MaxEE=2", aes(x = 0, y = log10(2), hjust = 0, vjust = 0), color = "red") +
      geom_text(label = "MaxEE=3", aes(x = 0, y = log10(3), hjust = 0, vjust = 0), color = "red") +
      geom_text(label = "MaxEE=5", aes(x = 0, y = log10(5), hjust = 0, vjust = 0), color = "red") +
      geom_text(label = "MaxEE=7", aes(x = 0, y = log10(7), hjust = 0, vjust = 0), color = "red") +
      labs(x = "Reads position", y = "EE = sum(10^(-Q/10)) log10") #+
    # coord_cartesian(ylim = c(log10(min(df.tmp[complete.cases(df.tmp), ]$A)), log10(max(df.tmp[complete.cases(df.tmp), ]$A))))
    return(q_r)
  }

  ### MAIN
  p_r1 <- qualPlot(df.tmp = fastqTmp(fastq.r1)) +
    ggtitle(name1) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  p_r2 <- qualPlot(df.tmp = fastqTmp(fastq.r2)) + ggtitle(name2) +
    ggtitle(name2) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  q_r1 <- maxEEplot(df.tmp = fastqTmp(fastq.r1))
  q_r2 <- maxEEplot(df.tmp = fastqTmp(fastq.r2))

  return(plot((p_r1 + p_r2) / (q_r1 + q_r2), cols = 2))
}




#' Summarise phyloseq table
#'
#' @param physeq
#' @param Rank
#' @param GroupBy
#'
#' @return
#'
#' @examples
summarize_taxa <- function(physeq, Rank, GroupBy = NULL) {
  Rank <- Rank[1]
  if (!Rank %in% rank_names(physeq)) {
    message(
      "The argument to `Rank` was:\n", Rank,
      "\nBut it was not found among taxonomic ranks:\n",
      paste0(rank_names(physeq), collapse = ", "), "\n",
      "Please check the list shown above and try again."
    )
  }
  if (!is.null(GroupBy)) {
    GroupBy <- GroupBy[1]
    if (!GroupBy %in% sample_variables(physeq)) {
      message(
        "The argument to `GroupBy` was:\n", GroupBy,
        "\nBut it was not found among sample variables:\n",
        paste0(sample_variables(physeq), collapse = ", "), "\n",
        "Please check the list shown above and try again."
      )
    }
  }
  # Start with fast melt
  mdt <- fast_melt(physeq)
  if (!is.null(GroupBy)) {
    # Add the variable indicated in `GroupBy`, if provided.
    sdt <- data.table(
      SampleID = sample_names(physeq),
      var1 = get_variable(physeq, GroupBy)
    )
    setnames(sdt, "var1", GroupBy)
    # Join
    setkey(sdt, SampleID)
    setkey(mdt, SampleID)
    mdt <- sdt[mdt]
  }
  # Summarize
  summarydt <- mdt[, list(totalRA = sum(RelativeAbundance)),
    by = c(Rank, GroupBy)
  ]
  return(summarydt)
}


# Create mismatch ---------------------------------------------------------


#' Create_mismatch
#'
#' @param dna
#' @param dist
#' @param ...
#'
#' @return
#'
#' @examples
create_mismatch <- function(dna, dist, ...) {
  all_bases <- c("A", "T", "C", "G")
  l <- tstrsplit(dna, "", fixed = TRUE)
  l <- lapply(l, function(x) all_bases)
  r <- Reduce(paste0, do.call(CJ, l))
  return(r[which(stringdist(dna, r, method = "hamming") <= dist)])
}



# Convert to proportions --------------------------------------------------



#' Convert phyloseq table to proportions
#'
#' @param x
#' @param thresh
#' @param na_rm
#' @param ...
#'
#' @return
#'
#' @examples
proportions <- function(x, thresh = NA, na_rm = FALSE, ...) {
  xprop <- (x / sum(x)) # Convert to proportions
  xprop[xprop <= thresh] <- NA ## remove taxa under this level
  xprop2 <- (xprop / sum(xprop, na.rm = na_rm))
  return(xprop2)
}


# ps_to_fasta --------------------------------------------------------------


#' outputs a FASTA file from a phyloseq object
#'
#' @description This function outputs a FASTA-formatted text file from a \code{phyloseq} object.
#' This code was modified from \code{reltools} package https://github.com/DanielSprockett/reltools by Daniel Sprocket
#'
#' @param ps A \code{phyloseq} object that contains \code{\link[phyloseq]{refseq}}.
#' If there the \code{refseq} slot is not filled, this function will try pull the
#' sequences from \code{\link[phyloseq]{get_taxa}}
#'
#' @param file (optional) A file name that ends in ".fasta" or ".fa".
#' If a file name is not supplied, the file will be named after the phyloseq object.
#'
#' @param rank (optional) A taxonomic rank from the \code{\link[phyloseq]{tax_table}} which will be used to name the sequences.
#' If no rank is supplied, samples will be named \code{ASV_#}
#'
#' @param width (Default 1000) The number of characters in each fasta line before wrapping occurs
#'
#' @param ... (Optional) Any further paramaters to be passed to \code{writeXStringSet}
#'
#' @return This function saves a FASTA-formatted text file from the input \code{phyloseq} object.
#' @export
#' @examples
#' save_fasta(ps)
#' save_fasta(ps = ps, file = "sequences.fasta", rank = "Genus")

ps_to_fasta <- function(ps = ps, out.file = NULL, rank = NULL, width = 1000, ...){

  if (is.null(ps)){
    message("Phyloseq object not found.")
  }

  if (is.null(out.file)){
    out.file <- paste0(deparse(substitute(ps)), ".fasta")
  }

  if (!is.null(refseq(ps, errorIfNULL = FALSE))){
    seqs <- DNAStringSet(as.vector(refseq(ps)))
  } else{
    message("refseq() not found. Using taxa names for sequences.")
    if (sum(grepl("[^ACTG]", rownames(tax_table(ps)))) > 0){
      stop("Error: Taxa do not appear to be DNA sequences.")
    }
    seqs <- DNAStringSet(colnames(get_taxa(ps)))
  }

  if (is.null(rank) || !rank %in% rank_names(ps)){
    message("Rank not found. Naming sequences sequentially (i.e. ASV_#).")
    names(seqs) <- paste0("ASV_", 1:ntaxa(ps))
  } else {
    names(seqs) <- make.unique(unname(tax_table(ps)[,rank]), sep = "_")
  }

  writeXStringSet(seqs, filepath = out.file, width=width, ... = ...)
  message(paste0(ntaxa(ps), " sequences written to <", out.file, ">."))
}


