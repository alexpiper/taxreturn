##Could turn this into a function with different inputs!

#Add a default model for clean_seqs

#First part - PHMM
clean_seqs <- function(x, model, minscore = 100, minamplen = 50,
                       maxamplen = 500,shave=FALSE, maxNs = 0.02, cores = 1,
                       quiet = FALSE){

  #Change to DNAbin
  x <- as.DNAbin(x)

  #Define PHMM function

  filt_phmm <- function(s, model, minscore, minamplen, maxamplen){
    s <- s[!s %in% as.raw(c(2, 4))]
    vit <- aphid::Viterbi(model, s, odds = TRUE, type = "semiglobal")

    if(vit$score < minscore) return(NULL)
    path <- vit$path
    match1 <- match(1, path)
    match2 <- match(1, rev(path))

    if(is.na(match1) | is.na(match2)) return(NULL)

    newlength <- length(vit$path) - match1 - match2 + 2
    if(newlength < minamplen | newlength > maxamplen) return(NULL)

    attr(s, "score") <- vit$score
    return(s)
  }

  nseq <- length(x)

  if(inherits(cores, "cluster")){
    x <- parallel::parLapply(cores, x, filt_phmm,  model, minscore, minamplen, maxamplen)
  }else if(cores == 1){
    x <- lapply(x, filt_phmm, model, minscore, minamplen, maxamplen)
  }else{
    navailcores <- parallel::detectCores()
    if(identical(cores, "autodetect")) cores <- navailcores - 1
    if(!(mode(cores) %in% c("numeric", "integer"))) stop("Invalid 'cores'")
    # if(cores > navailcores) stop("Number of cores is more than available")
    if(cores > 1){
      if(!quiet) cat("Multithreading with", cores, "cores\n")

      cores <- parallel::makeCluster(cores, outfile="out.txt")
      junk <- clusterEvalQ(cores, sapply(c("bold","taxizedb","tidyverse","rentrez","Biostrings","biofiles"), require, character.only = TRUE)) #Discard result

      x <- parallel::parLapply(cores, filt_phmm,  model, minscore, minamplen, maxamplen)
      parallel::stopCluster(cores)
    }else{
      x <- lapply(x, filt_phmm, model, minscore, minamplen, maxamplen)
    }
  }

  discards <- sapply(x, is.null)
  nseq <- sum(!discards)


  if(nseq > 0){
    if(!quiet) cat("Retained", nseq, "sequences after alignment to PHMM\n")
    scores <- unlist(lapply(x, function(s) attr(s, "score")), use.names = FALSE)
    x <- x[!discards]
    x <- as.DNAbin(as.character.DNAbin(x))

  }else{
    if(!quiet) cat("None of the sequences met PHMM specificity criteria. Returning NULL\n")
    x <- NULL
  }
  if(!quiet) cat("Filtering ambiguous sequences\n")
  discards <- sapply(x, function(s) sum(s == 0xf0)/length(s)) > maxNs
  x <- subset.DNAbin(x, subset = !discards)
  if(!quiet) cat(length(x), "sequences retained after applying ambiguity filter\n")
  if(shave==TRUE){
   x <- shave(x, model, direction = "both", cores)
   if(!quiet) cat("Overhanging sequence data shaved from PHMM\n")
  }
  if(!quiet) cat("Done\n")
  return(x)
}

## Prune groups function - See if there is a way to removed by length - Can we sort the groups by seq length then start from bottom
#Add a discardby=Random, or discardby=Length

prune_groups <- function(x, maxGroupSize=5, quiet=FALSE){
  groups <- names(x) %>%
    str_split_fixed(";",n=2) %>%
    as_tibble() %>%
    separate(V1,into=c("acc","taxid"))%>%
    pull(taxid)
  groupCounts <- table(groups)
  u_groups <- names(groupCounts)

  remove <- logical(length(x))
  for (i in which(groupCounts > maxGroupSize)) {
    index <- which(groups == u_groups[i])
    keep <- sample(
      length(index),
      maxGroupSize
    )
    remove[index[-keep]] <- TRUE
  }
  x <- x[!remove]
  if(!quiet) cat(paste0(sum(remove), " sequences pruned from over-represented groups"))
  return(x)
}

#filter_taxa <- function(file,minlength,maxlength,unique=TRUE,binomials=TRUE, removeterms){


###Filter any further erroneous or insufficiently identified sequences
#
#bold_names <- as_tibble(unlist(getAnnot(bold)))
#bold_names <- bold_names %>%
#  dplyr::filter(!str_detect(value, fixed("sp."))) %>%
#  dplyr::filter(!str_detect(value, fixed("aff."))) %>%
#  dplyr::filter(!str_detect(value, fixed("nr."))) %>%
#  dplyr::filter(!str_detect(value, fixed("cf."))) %>%
#  dplyr::filter(!str_detect(value, fixed("nom."))) %>%
#  dplyr::filter(!str_detect(value, fixed("nud."))) %>%
#  dplyr::filter(!str_detect(value, fixed("environment"))) %>%
#  dplyr::filter(!str_detect(value, fixed("undescribed"))) %>%
#  dplyr::filter(!str_detect(value, fixed("unverified"))) %>%
#  dplyr::filter(!str_detect(value, fixed("uncultured"))) %>%
#  dplyr::filter(!str_detect(value, fixed("unidentif"))) %>%
#  dplyr::filter(!str_detect(value, fixed("Bacterium"))) %>%
#  dplyr::filter(!str_detect(value, fixed("wolbachia"))) %>%
#  dplyr::filter(!str_detect(value, fixed("symbiont"))) %>%
#  dplyr::filter(!str_detect(value, fixed("Bacterium"))) %>%
#  dplyr::filter(!str_detect(value, fixed("NA"))) %>%
#  dplyr::filter(!str_detect(value, fixed("error")))
#
#rm_keywords <- bold_names$value
#
#bold <- bold[getAnnot(bold) %in% rm_keywords]
#bold_filtered$keywords <- length(bold)
#
##Write out filtered fasta
#bold_names <- getName(bold)
#write.fasta(bold, bold_names, paste0("bold_tempfilt1.fa"), as.string=FALSE, nbchar=100)
#
##test
#bold <- readFASTA("bold_tempfilt1.fa")
#
##trim to amplicon - BF1 - BR1
##this is probably better to replace with a PHMM directly for the sequencec
#amplicon <- virtualPCR(bold, up = "ACWGGWTGRACWGTNTAYCC",down= "ARYATDGTRATDGCHCCDGC",cores=3, rcdown = TRUE, trimprimers = TRUE)
#writeFASTA(amplicon,"bold_trimmed.fa")
#
#bold_filtered$trimmed <- length(amplicon)
#
#rm(amplicon)
#rm(bold)
#
