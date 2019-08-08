
#' Clean Sequences with PHMM
#'
#' @param x
#' @param model
#' @param minscore
#' @param shave
#' @param maxNs
#' @param cores
#' @param quiet
#'
#'
#' @import tidyverse
#' @import aphid
#' @import insect
#' @import Biostrings
#' @import ape
#' @import stringr
#' @import parallel
#' @import pbapply
#'
#' @return
#' @export
#'
#' @examples
clean_seqs <- function(x, model, minscore = 100, shave=TRUE, maxNs = 0.02, cores = 1,
                       quiet = FALSE, progress=FALSE){
  time <- Sys.time() # get time
  #Convert to DNAbin
  if(!is(x,"DNAbin")){ x <- ape::as.DNAbin(x)}

  #Define PHMM function

  filt_phmm <- function(s, model, minscore, minamplen, maxamplen){

    #.packages=c("aphid","insect","ape")
    s <- s[!s %in% as.raw(c(2, 4))]
    vit <- aphid::Viterbi(model, s, odds = TRUE, type = "semiglobal",cpp=TRUE,residues="DNA")

    if(vit$score < minscore) return(NULL)
    path <- vit$path
    match1 <- match(1, path)
    match2 <- match(1, rev(path))

    if(is.na(match1) | is.na(match2)) return(NULL)

    if(shave==TRUE){
      ntoshavef <- match(c(0, 1), rev(vit$path)) - 1
      ntoshavef <- min(ntoshavef[!is.na(ntoshavef)])
      last <- length(s) - ntoshavef
      begin <- match(c(0, 1), vit$path)
      begin <- min(begin[!is.na(begin)])
      s <- s[begin:last]
    }
    attr(s, "score") <- vit$score
    return(s)
  }

  nseq <- length(x)

  if(cores == 1 && progress == TRUE){
    x <- pbapply::pblapply(x, filt_phmm, model, minscore)
    } else if(cores == 1 && progress == FALSE){
      x <- lapply(x, filt_phmm, model, minscore)
  } else if(cores > 1 && progress == TRUE){
    stop("Progress bar currently not supported for multithreading")
  }  else{
    navailcores <- parallel::detectCores()
    if(identical(cores, "autodetect")) cores <- navailcores - 1
    if(!(mode(cores) %in% c("numeric", "integer"))) stop("Invalid 'cores'")
     if(cores > navailcores) stop("Number of cores is more than available")

    if(cores > 1){
      if(!quiet) cat("Multithreading with", cores, "cores\n")

      cores <- parallel::makeCluster(cores, outfile="out.txt")
      #parallel::clusterExport(cores, c("model", "minscore"))
      junk <- parallel::clusterEvalQ(cores, sapply(c("aphid","insect","ape"), require, character.only = TRUE)) #Discard result

      x <- parallel::parLapply(cores,x, filt_phmm,  model, minscore)
      parallel::stopCluster(cores)
    }else{
      x <- lapply(x, filt_phmm, model, minscore)
    }
  }

  discards <- sapply(x, is.null)
  nseq <- sum(!discards)

  if(nseq > 0){
    if(!quiet) cat("Retained", nseq, "sequences after alignment to PHMM\n")
    scores <- unlist(lapply(x, function(s) attr(s, "score")), use.names = FALSE)
    x <- x[!discards]
    x <- ape::as.DNAbin(ape::as.character.DNAbin(x))

  }else{
    if(!quiet) cat("None of the sequences met PHMM specificity criteria. Returning NULL\n")
    x <- NULL
  }
  if(!quiet) cat("Filtering ambiguous sequences\n")
  discards <- sapply(x, function(s) sum(s == 0xf0)/length(s)) > maxNs
  x <- insect::subset.DNAbin(x, subset = !discards)
  if(!quiet) cat(length(x), "sequences retained after applying ambiguity filter\n")
  if(!quiet) cat("Bases overhanging PHMM shaved from alignment\n")
  if(!quiet) cat("Done\n")
  time <- Sys.time() - time
  if (!quiet) (message(paste0("finished in ", format(time, digits=2))))
  return(x)
}


# Prune group sizes -------------------------------------------------------


## Prune groups function - See if there is a way to removed by length - Can we sort the groups by seq length then start from bottom
#Add a discardby=Random, or discardby=Length

#' Prune group sizes
#'
#' @param x
#' @param maxGroupSize
#' @param quiet
#'
#' @return
#' @export
#'
#' @import tidyverse
#' @examples
prune_groups <- function(x, maxGroupSize=5,removeby="random", quiet=FALSE){
  groups <- names(x) %>%
    str_split_fixed(";",n=2) %>%
    as_tibble() %>%
    separate(V1,into=c("acc","taxid"))%>%
    pull(taxid)
  groupCounts <- table(groups) # Count number of seqs per group
  u_groups <- names(groupCounts) #Get unique groups

  remove <- logical(length(x))
  if(removeby=="random"){
    for (i in which(groupCounts > maxGroupSize)) {
      index <- which(groups == u_groups[i])
      keep <- sample( # Take random sample
        length(index),
        maxGroupSize
      )
      remove[index[-keep]] <- TRUE
    }
  }else if (removeby=="length"){
    for (i in which(groupCounts > maxGroupSize)) {

      index <- which(groups == u_groups[i])

      rem <- lengths(x[index])
      names(rem) <- index
      rem <- sort(rem, decreasing = TRUE)

      keep <- as.integer(names(rem[1:maxGroupSize]))
      remove[index[!index %in% keep]] <- TRUE

    }
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
