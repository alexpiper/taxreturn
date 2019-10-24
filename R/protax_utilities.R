
# Define reference sequences ----------------------------------------------


protax_get_reference_sequences <- function(install=NULL, taxonomy, trainseq2taxid, output, method="random", levels=NULL, max_per_level=NULL, seed=1){

    install <- paste0(install, "/scripts/")
    taxonomy <- normalizePath(paste0(getwd(),taxonomy))
    trainseq2taxid <- normalizePath(paste0(getwd(),trainseq2taxid))
    output <- paste0(getwd(),output)

    if (method == "random"){
      script <- paste0(install, "/get_reference_sequences_random.pl")
    } else if (method == "all"){
      script <- paste0(install, "/get_reference_sequences_all.pl")
      stop("Warning: only random selection of reference sequences is currently working")
    } else if (method == "clustering"){
      script <- paste0(install, "/get_reference_sequences_clustering.pl")
      stop("Warning: only random selection of reference sequences is currently working")
    }

    #Need to define how many levels contained within the reference taxonomy tree
    if(is.null(levels)){stop("Number of levels in input taxonomy needs to be defined")}

    #max = how many sequences are allowed for each level in the input taxonomy
    if(length(max_per_level == 1)){
      max_per_level = paste(rep(max_per_level, levels), collapse = " ")
    }

    args <- paste(normalizePath(script) ,taxonomy, levels, max_per_level, trainseq2taxid, seed, output)

    print(paste0(args," /n"))
    #Run perl script
    system2("perl", args = args)

}


# Generate training data --------------------------------------------------


protax_generate_training_data <- function(install=NULL, taxonomy, trainseq2taxid, repseqs, output, samples=NULL, seed=1, include_ignored=FALSE){

  install <- paste0(install, "/scripts/")
  taxonomy <- normalizePath(paste0(getwd(),taxonomy))
  trainseq2taxid <- normalizePath(paste0(getwd(),trainseq2taxid))
  repseqs <- normalizePath(paste0(getwd(),repseqs))
  output <- paste0(getwd(),output)

  script <- paste0(install, "/generate_training_data2.pl")

  if(is.null(samples)){stop("Number of samples to be drawn form input taxonomy needs to be defined")}

  #option to ignore certain sequences when calculating probabilities
  if(include_ignored == FALSE){
    ignore="no"
    }else (stop("Including sequences to be ignored is not currently implemented"))

  args <- paste(normalizePath(script), taxonomy, trainseq2taxid, repseqs, samples, seed, ignore, output)

  print(paste0(args," /n"))
  #Run perl script
  system2("perl", args = args)

}



# Create X matrix ---------------------------------------------------------

#create_xdata_seqsimfile.pl my_trainsamples.txt example_taxonomy.txt example_trainseqid2taxname.txt my_rseqs.txt example_trainseqsim.txt my_trainxdata.txt

protax_create_xmatrix <- function(install=NULL, repsamples, taxonomy, trainseq2taxid, repseqs, seqsim, output){

  install <- paste0(install, "/scripts/")
  taxonomy <- normalizePath(paste0(getwd(), taxonomy))
  repsamples <- normalizePath(paste0(getwd(), repsamples))
  trainseq2taxid <- normalizePath(paste0(getwd(), trainseq2taxid))
  repseqs <- normalizePath(paste0(getwd(), repseqs))
  seqsim <- normalizePath(paste0(getwd(), seqsim))
  output <- paste0(getwd(), output)

  script <- paste0(install, "/create_xdata_seqsimfile.pl")

  args <- paste(normalizePath(script), repsamples, taxonomy, trainseq2taxid, repseqs, seqsim, output)

  print(paste0(args," /n"))
  #Run perl script
  system2("perl", args = args)

}


# Classify sequences ------------------------------------------------------

#This should be rewritten to take in a fasta and output the seqsims using kmer, and the seqids for classification

protax_classify <- function(install=NULL, seqids, seqsim, taxonomy, trainseq2taxid,
                            repseqs, mcmc, output, type="map", n_outcomes=1,
                            threshold =0.01, add_tax=TRUE, validation=FALSE, verbose=TRUE,
                            nodeprob=FALSE){

  if (!type=="map"){stop("Including sequences to be ignored is not currently implemented")}
  if (validation==FALSE){
    val = 0
  }else if (validation ==TRUE){val = 1}

  if (verbose==FALSE){
    verbose = 0
  }else if (verbose ==TRUE){verbose = 1}

  install <- paste0(install, "/scripts/")
  seqids <- normalizePath(paste0(getwd(), seqids))
  seqsim <- normalizePath(paste0(getwd(), seqsim))
  taxonomy <- normalizePath(paste0(getwd(), taxonomy))
  trainseq2taxid <- normalizePath(paste0(getwd(), trainseq2taxid))
  repseqs <- normalizePath(paste0(getwd(), repseqs))
  mcmc <- normalizePath(paste0(getwd(), mcmc))
  output <- paste0(getwd(), output)

  if(nodeprob==FALSE){
  script <- paste0(install, "/classify_seqsimfile.pl")
  } else if (nodeprob==TRUE){script <- paste0(install, "/nodeprob_seqsimfile.pl")}

  args <- paste(normalizePath(script), seqids, taxonomy, trainseq2taxid, repseqs, mcmc, type, seqsim, n_outcomes, threshold, output, val, verbose)

  print(paste0(args," /n"))
  #Run perl script
  system2("perl", args = args) # Direct classification back as an R object

 #Add tax - could probably be done manually in R
  if(add_tax==TRUE){
    script2 <- paste0(install, "/add_taxonomy_info.pl")

    args2 <- paste(normalizePath(script2), taxonomy, output)

    print(paste0(args2," /n"))

    x <- read_tsv(system2("perl", args = args2, stdout=TRUE), col_names=FALSE) %>%
      magrittr::set_colnames(c("query", "nodeid", "prob","level", "nodename"))
    return(x)
  }

}



