
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
clean_seqs <- function(x, model, minscore = 100, shave = TRUE, maxNs = 0, cores = 1,
                       quiet = FALSE, progress = FALSE, ...) {
  time <- Sys.time() # get time

  # Convert to DNAbin
  if (!is(x, "DNAbin")) {
    x <- ape::as.DNAbin(x)
    if (all(is.na(ape::base.freq(x)))) {stop("Error: Object is not coercible to DNAbin \n")}
  }
  if (!is(model, "PHMM")) {
    stop("Model needs to be a PHMM object")
  }
  # Define PHMM function

  filt_phmm <- function(s, model, minscore, minamplen, maxamplen, ...) {

    s <- s[!s %in% as.raw(c(2, 4))]
    vit <- aphid::Viterbi(model, s, odds = TRUE, type = "semiglobal", cpp = TRUE, residues = "DNA")

    if (vit$score < minscore) {
      return(NULL)
    }
    path <- vit$path
    match1 <- match(1, path)
    match2 <- match(1, rev(path))

    if (is.na(match1) | is.na(match2)) {
      return(NULL)
    }

    if (shave == TRUE) {
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

  if (cores == 1 && progress == TRUE) {
    x <- pbapply::pblapply(x, filt_phmm, model, minscore)
  } else if (cores == 1 && progress == FALSE) {
    x <- lapply(x, filt_phmm, model, minscore)
  } else if (cores > 1 && progress == TRUE) {
    stop("Progress bar currently not supported for multithreading")
  } else {
    navailcores <- parallel::detectCores()
    if (identical(cores, "autodetect")) cores <- navailcores - 1
    if (!(mode(cores) %in% c("numeric", "integer"))) stop("Invalid 'cores'")
    if (cores > navailcores) stop("Number of cores is more than available")

    if (cores > 1) {
      if (!quiet) cat("Multithreading with", cores, "cores\n")

      cores <- parallel::makeCluster(cores, outfile = "out.txt")
      # parallel::clusterExport(cores, c("model", "minscore"))
      junk <- parallel::clusterEvalQ(cores, sapply(c("aphid", "insect", "ape"), require, character.only = TRUE)) # Discard result

      x <- parallel::parLapply(cores, x, filt_phmm, model, minscore)
      parallel::stopCluster(cores)
    } else {
      x <- lapply(x, filt_phmm, model, minscore)
    }
  }

  discards <- sapply(x, is.null)
  nseq <- sum(!discards)

  if (nseq > 0) {
    if (!quiet) cat("Retained", nseq, "sequences after alignment to PHMM\n")
    scores <- unlist(lapply(x, function(s) attr(s, "score")), use.names = FALSE)
    x <- x[!discards]
    x <- ape::as.DNAbin(ape::as.character.DNAbin(x))
  } else {
    if (!quiet) cat("None of the sequences met PHMM specificity criteria. Returning NULL\n")
    x <- NULL
  }
  if (!quiet) cat("Filtering ambiguous sequences\n")
  discards <- sapply(x, function(s) sum(s == 0xf0) / length(s)) > maxNs
  x <- insect::subset.DNAbin(x, subset = !discards)
  if (!quiet) cat(length(x), "sequences retained after applying ambiguity filter\n")
  if (!quiet) cat("Bases overhanging PHMM shaved from alignment\n")
  if (!quiet) cat("Done\n")
  time <- Sys.time() - time
  if (!quiet) (message(paste0("finished in ", format(time, digits = 2))))
  return(x)
}


# Prune group sizes -------------------------------------------------------

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
prune_groups <- function(x, maxGroupSize = 5, dedup = TRUE, discardby = "random", quiet = FALSE) {
  if (dedup) {
    dup <- length(x)
    x <- insect::subset.DNAbin(x, subset = !insect::duplicated.DNAbin(x, point = TRUE))
    if (!quiet) cat(paste0((dup - length(x)), " duplicate sequences removed \n"))
  }

  groups <- names(x) %>%
    str_split_fixed(";", n = 2) %>%
    as_tibble() %>%
    separate(V1, into = c("acc", "taxid")) %>%
    pull(taxid)
  groupCounts <- table(groups) # Count number of seqs per group
  u_groups <- names(groupCounts) # Get unique groups

  remove <- logical(length(x))
  if (discardby == "random") {
    for (i in which(groupCounts > maxGroupSize)) {
      index <- which(groups == u_groups[i])
      keep <- sample( # Take random sample
        length(index),
        maxGroupSize
      )
      remove[index[-keep]] <- TRUE
    }
  } else if (discardby == "length") {
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
  if (!quiet) cat(paste0(sum(remove), " sequences pruned from over-represented groups"))
  return(x)
}


# Resolve synonyms -------------------------------------------------------

#' Resolve taxonomic synonyms
#' @description This function takes a DNAbin object, or a list of species and uses the Global Biodiversity Information Facility (GBIF) to resolve taxonomic synonyms
#'
#' @param x
#' @param subspecies
#' @param quiet
#' @param missing
#' @param higherrank
#' @param fuzzy
#'
#' @return
#' @export
#'
#' @import tidyverse
#' @examples
resolve_taxonomy <- function(x, subspecies = FALSE, quiet = TRUE, missing = "ignore", higherrank = FALSE, fuzzy = TRUE) {
  time <- Sys.time() # get time

  if (quiet == TRUE) {
    verbose <- FALSE
  } else {
    (verbose <- TRUE)
  }

  # Check type of input
  if (is(x, "DNAbin")) {
    message("Input is DNAbin")
    is.seq <- TRUE
  } else  if (is(x, "DNAStringSet")| is(x, "DNAString")) {
    message("Input is DNAStringSet, converting to DNAbin")
    x <- ape::as.DNAbin(x)
    is.seq <- TRUE
  } else  if (is(x, "character")) {
    message("Input is character vector, resolving Genus species binomials")
    is.seq <- FALSE
  }

  # if input has sequences, get names
  if (is.seq == TRUE) {
  query <- names(x) %>%
    stringr::str_split_fixed(";", n = 2) %>%
    tibble::as_tibble() %>%
    tidyr::separate(col = V1, into = c("acc", "taxid"), sep = "\\|") %>%
    dplyr::rename(query = V2)
  } else if (is.seq ==FALSE ) {
    query <- data.frame(query= x)
    query$taxid <- "NA" # Add dummy columns
    query$acc <- "NA"
    }

  if (verbose == TRUE) {
    message(paste0(
      "resolving synonyms for ", length(unique(query$query)),
      " Unique taxa, estimated time: ", signif((length(unique(query$query)) * 0.5) / 3600, digits = 3), " hours" #need better estimator
    ))
  }
  out <- resolve_gbif(unique(query$query), subspecies = subspecies, verbose = verbose, higherrank = higherrank, fuzzy = fuzzy, resolve_taxonomy = TRUE)

  out <- out %>%
    dplyr::filter(synonym == TRUE) %>%
    dplyr::filter(!scientificName == "Curculio bimaculatus") %>% # This taxa is breaking function
    dplyr::mutate(taxidnew = taxizedb::name2taxid(scientificNameStd)) %>%
    dplyr::mutate(query = as.character(scientificName)) %>%
    dplyr::mutate(scientificNameStd = as.character(scientificNameStd)) %>%
    dplyr::mutate_all(as.character())

  if(nrow(out) > 0) {
    if (missing == "ignore") { # In cases where the updated synonym does not have a taxonomic ID in the NCBI database, update name but keep old taxid

      query <- query %>%
        dplyr::left_join(out, by = "query") %>%
        dplyr::mutate(query = case_when(
          !is.na(scientificNameStd) ~ scientificNameStd,
          is.na(scientificNameStd) ~ query
        )) %>% # Catch all for anything that is not above case
        dplyr::mutate(taxid = case_when(
          !is.na(taxidnew) ~ taxidnew,
          TRUE ~ taxid
        )) %>%
        dplyr::select(acc, taxid, query)
    } else if (missing == "keepold") { # In cases where the updated synonym does not have a taxonomic ID in the NCBI database, keep old column

      query <- query %>%
        dplyr::left_join(out, by = "query") %>%
        dplyr::mutate(query = case_when(
          !is.na(scientificNameStd) & !is.na(taxidnew) ~ scientificNameStd,
          TRUE ~ query
        )) %>% # Catch all for anything that is not above case
        dplyr::mutate(taxid = case_when(
          !is.na(scientificNameStd) & !is.na(taxidnew) ~ taxidnew,
          TRUE ~ taxid
        )) %>%
        dplyr::select(acc, taxid, query)
    } else if (missing == "remove") { # In cases where the updated synonym does not have a taxonomic ID in the NCBI database, remvoe

      query <- query %>%
        dplyr::left_join(out, by = "query") %>%
        dplyr::mutate(keepcol = case_when(
          !is.na(scientificNameStd) & is.na(taxidnew) ~ FALSE,
          !is.na(scientificNameStd) & !is.na(taxidnew) ~ TRUE,
          is.na(scientificNameStd) & !is.na(query) ~ TRUE,
          TRUE ~ TRUE
        )) %>% # Catch all for anything that is not above case
        dplyr::filter(keepcol == TRUE) %>%
        dplyr::mutate(query = case_when(
          !is.na(scientificNameStd) ~ scientificNameStd,
          TRUE ~ query
        )) %>% # Catch all for anything that is not above case
        dplyr::mutate(taxid = case_when(
          !is.na(taxidnew) ~ taxidnew,
          TRUE ~ taxid
        )) %>%
        dplyr::select(acc, taxid, query)
    }

    if(is.seq == TRUE) {
    query <- query %>%
      tidyr::unite(col = V1, c("acc", "taxid"), sep = "|")

      names(x) <- paste(query$V1, query$query, sep = ";")
    } else if (is.seq == FALSE) {
        x <- query$query
      }

  } else (message("No synonyms detected"))


  time <- Sys.time() - time
  if (!quiet) {message(paste0("resolved ", nrow(out), " synonyms in ", format(time, digits = 2)))}
  return(x)
}



### Get GBIF Taxonomy

# modified from traitdataform pacakge https://github.com/EcologicalTraitData/traitdataform/issues/35

# COuld remove subspecies
# COuld te

resolve_gbif <- function(x, subspecies = TRUE, higherrank = TRUE, verbose = FALSE,
                                     fuzzy = TRUE, conf_threshold = 90, resolve_taxonomy = TRUE) {
  matchtype <- status <- confidence <- NULL

  ## Get GBIF data - this needs to be sped up majorly, would be nice to move to taxize::db
  # temp <- taxize::get_gbifid_(x, messages = verbose)

  temp <- x %>%
    purrr::map(safely(taxize::get_gbifid_)) %>%
    purrr::map("result") %>%
    flatten()

  i <- 1
  for (i in 1:length(temp)) {
    warning_i <- ""
    synonym_i <- FALSE
    if (nrow(temp[[i]]) == 0) {
      warning_i <- paste("No matching species concept!")
      temp[[i]] <- data.frame(
        scientificName = x[i], matchtype = "NONE",
        status = "NA", rank = "species"
      )
    }
    if (!fuzzy & nrow(temp[[i]]) > 0) {
      temp[[i]] <- subset(temp[[i]], matchtype != "FUZZY")
      if (nrow(temp[[i]]) == 0) {
        warning_i <- paste(warning_i, "Fuzzy matching might yield results.")
      }
    }
    if (!is.null(conf_threshold) & nrow(temp[[i]]) > 0) {
      temp[[i]] <- subset(temp[[i]], confidence >= conf_threshold)
      if (nrow(temp[[i]]) == 0) {
        temp[[i]] <- data.frame(
          scientificName = x[i],
          matchtype = "NONE", status = "NA", rank = "species"
        )
        warning_i <- paste(warning_i, "No match! Check spelling or lower confidence threshold!")
      }
    }
    if (any(temp[[i]]$status == "ACCEPTED")) {
      temp[[i]] <- subset(temp[[i]], status == "ACCEPTED")
      temp[[i]] <- subset(temp[[i]], temp[[i]]$confidence ==
        max(temp[[i]]$confidence))
      if (nrow(temp[[i]]) > 1) {
        temp[[i]] <- temp[[i]][1, ]
        warning_i <- paste(warning_i, "Selected first of multiple equally ranked concepts!")
      }
    }
    if (!any(temp[[i]]$status == "ACCEPTED") & any(temp[[i]]$status ==
      "SYNONYM")) {
      if (resolve_taxonomy) {
        keep <- temp[i]
        if (!is.null(temp[[i]]$species) && !is.na(temp[[i]]$species)) {
          temp[i] <- taxize::get_gbifid_(temp[[i]]$species[which.max(temp[[i]]$confidence)],
            messages = verbose
          )
        } else if (is.null(temp[[i]]$species)) {
          newspp <- str_split_fixed(names(temp[i]), pattern = " ", n = 2)
          temp[[i]]$species <- paste(temp[[i]]$genus, newspp[1, 2])
          temp[i] <- taxize::get_gbifid_(temp[[i]]$species[which.max(temp[[i]]$confidence)],
            messages = verbose
          )
          if (nrow(temp[[i]]) < 1) {
            temp[i] <- keep
          }
        }

        if (temp[[i]][1, ]$status == "ACCEPTED" & !temp[[i]][1, ]$matchtype == "HIGHERRANK") {
          # if (temp[[i]][1, ]$status == "ACCEPTED") {
          temp[[i]] <- subset(temp[[i]], status == "ACCEPTED") # , matchtype == "EXACT" &
          temp[[i]] <- subset(temp[[i]], temp[[i]]$confidence ==
            max(temp[[i]]$confidence))
          if (nrow(temp[[i]]) > 1) {
            temp[[i]] <- temp[[i]][1, ]
            warning_i <- paste(warning_i, "Selected first of multiple equally ranked concepts!")
          }
          warning_i <- paste(warning_i, "A synonym was mapped to the accepted species concept!",
            sep = " "
          )
          synonym_i <- TRUE
        } else {
          status <- temp[[i]][1, ]$status
          temp[i] <- keep # Putting it back to before the call?
          if (nrow(temp[[i]]) > 1) {
            temp[[i]] <- temp[[i]][1, ]
            warning_i <- paste(warning_i, "Selected first of multiple equally ranked concepts!")
          }
          warning_i <- paste0(
            warning_i, " Resolved synonym '",
            temp[[i]]$species, "' is labelled '", status,
            "'. Clarification required!"
          )
        }
      } else {
        temp[[i]] <- subset(temp[[i]], status == "SYNONYM")
        temp[[i]] <- subset(temp[[i]], temp[[i]]$confidence ==
          max(temp[[i]]$confidence))
        warning_i <- paste(warning_i, "The provided taxon seems to be a synonym of '",
          temp[[i]]$species, "'!",
          sep = ""
        )
      }
    }
    if (all(temp[[i]]$status == "DOUBTFUL")) {
      temp[[i]] <- subset(temp[[i]], status == "DOUBTFUL")
      warning_i <- paste(warning_i, "Mapped concept is labelled 'DOUBTFUL'!")
      temp[[i]] <- subset(temp[[i]], temp[[i]]$confidence ==
        max(temp[[i]]$confidence))
      if (nrow(temp[[i]]) > 1) {
        temp[[i]] <- temp[[i]][1, ]
        warning_i <- paste(warning_i, "Selected first of multiple equally ranked concepts!")
      }
    }
    rankorder <- c(
      "kingdom", "phylum", "order", "class",
      "family", "genus", "species", "subspecies"
    )
    if (match(temp[[i]]$rank, rankorder) > 7 & !subspecies) {
      if (length(strsplit(
        as.character(temp[[i]]$canonicalname),
        " "
      )[[1]]) > 2) {
        temp[i] <- taxize::get_gbifid_(paste(strsplit(
          names(temp[i]),
          " "
        )[[1]][1:2], collapse = " "), messages = verbose)
        temp[[i]] <- subset(temp[[i]], temp[[i]]$confidence ==
          max(temp[[i]]$confidence))
        warning_i <- paste(warning_i, "Subspecies has been remapped to species concept!",
          sep = " "
        )
      } else {
        temp[[i]] <- data.frame(
          scientificName = x[i],
          matchtype = "NONE", rank = "subspecies"
        )
        warning_i <- paste(warning_i, "No mapping of subspecies name to species was possible!",
          sep = " "
        )
      }
    }
    if (temp[[i]]$matchtype == "HIGHERRANK") {
      if (higherrank) {
        temp[[i]] <- subset(temp[[i]], temp[[i]]$confidence ==
          max(temp[[i]]$confidence))
        warning_i <- paste(warning_i, "No matching species concept! Entry has been mapped to higher taxonomic level.")
      }
      else {
        temp[[i]] <- data.frame(
          scientificName = x[i],
          matchtype = "NONE", rank = "highertaxon"
        )
        warning_i <- paste(
          "No matching species concept!",
          warning_i
        )
      }
    }
    if (temp[[i]]$matchtype != "NONE") {
      temp[[i]] <- data.frame(
        scientificName = x[i], synonym = synonym_i,
        scientificNameStd = temp[[i]]$canonicalname,
        author = sub(paste0(
          temp[[i]]$canonicalname,
          " "
        ), "", temp[[i]]$scientificname), taxonRank = temp[[i]]$rank,
        confidence = temp[[i]]$confidence, kingdom = if (is.null(temp[[i]]$kingdom)) {
          NA
        } else {
          temp[[i]]$kingdom
        }, phylum = if (is.null(temp[[i]]$phylum)) {
          NA
        } else {
          temp[[i]]$phylum
        }, class = if (is.null(temp[[i]]$class)) {
          NA
        } else {
          temp[[i]]$class
        }, order = if (is.null(temp[[i]]$order)) {
          NA
        } else {
          temp[[i]]$order
        }, family = if (is.null(temp[[i]]$family)) {
          NA
        } else {
          temp[[i]]$family
        }, genus = if (is.null(temp[[i]]$genus)) {
          NA
        } else {
          temp[[i]]$genus
        }, taxonomy = "GBIF Backbone Taxonomy",
        taxonID = paste0(
          "http://www.gbif.org/species/",
          temp[[i]]$usagekey, ""
        ), warnings = NA
      )
    } else {
      temp[[i]] <- data.frame(scientificName = x[i], warnings = NA) # FAILING HERE?
    }
    temp[[i]]$warnings <- warning_i
    if (verbose & nchar(warning_i) >= 1) {
      warning(warning_i)
    }
  }
  out <- data.table::rbindlist(temp, fill = TRUE)
  class(out) <- c("data.frame", "taxonomy")
  return(out)
}


# Get Reading frame ---------------------------------------------------------------

#' Get Reading frame of sequences
#'
#' @param x Sequences in DNAStringset or DNAbin format
#' @param genetic.code A genetic code for the Amino acid translation. See all known codes at GENETIC_CODE_TABLE
#' Default is the invertebrate mitochondrial code 'SGC4'
#' @param forward
#' @param reverse
#' @param resolve_draws How draws should be resolved when multiple possible frames produce sequences with no stop codons.
#'
#' @return
#' @export
#'
#' @examples
get_reading_frame <- function(x, genetic.code = "SGC4", forward=TRUE, reverse=FALSE, resolve_draws="majority") {
  # Convert to DNAbin
  if (is(x, "DNAbin")) {
    x <- x %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
  }
  if(forward==TRUE) {
    F_frames <- lapply(1:3, function(pos) subseq(x, start=pos))
  }
  if(reverse==TRUE) {
    R_frames <- lapply(1:3, function(pos) subseq(reverseComplement(x), start=pos))
  }
  #Translate all reading frames
  suppressWarnings(translated <- lapply(F_frames, translate, genetic.code = getGeneticCode(genetic.code)))
  #select the reading frames that contain 0 stop codons, or return NA
  reading_frame <- vector("integer", length=length(x))
  for (i in 1:length(x)){
    fvec = c(str_count(as.character(translated[[1]][i]), "\\*"),
             str_count(as.character(translated[[2]][i]), "\\*"),
             str_count(as.character(translated[[3]][i]), "\\*"))
    if(sum(fvec==0)==1){
      reading_frame[i] <- which(fvec==0)
    } else if(sum(fvec==0)>1) {
      reading_frame[i] <- 0
    }else if(sum(fvec==0)==0) {
      reading_frame[i] <- NA
    }
  }
  if (resolve_draws == "majority") {
    reading_frame[reading_frame==0] <- reading_frame[which.max(tabulate(reading_frame))]
  } else if (resolve_draws == "remove") {
    reading_frame[reading_frame==0] <- NA
  }
  return(reading_frame)
}


# Codon_filter ------------------------------------------------------------

#' Filter sequences containing stop codons
#'
#' @param x Sequences in DNAStringset or DNAbin format
#' @param genetic.code A genetic code for the Amino acid translation. See all known codes at GENETIC_CODE_TABLE
#' Default is the invertebrate mitochondrial code 'SGC4'
#' @param forward
#' @param reverse
#'
#' @return
#' @export
#'
#' @examples
codon_filter <- function(x, genetic.code = "SGC4", forward=TRUE, reverse=FALSE){
  #Get reading frames
  frames <- get_reading_frame(x, genetic.code = genetic.code, forward = forward, reverse = reverse)

  out <- x[!is.na(frames)]
  message(paste0(length(x) - length(out), " Sequences containing stop codons removed"))
  return(out)
}



# Codon entropy  -----------------------------------------------------------

#' Codon entropy
#'
#' @param x
#' @param genetic.code
#' @param forward
#' @param reverse
#' @param codon.filter
#'
#' @return
#' @export
#'
#' @examples
codon_entropy <- function(x, genetic.code = "SGC4", forward=TRUE, reverse=FALSE, codon.filter = TRUE, method="ML") {
  if (is(x, "DNAbin")) {
    x <- x %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
  }
  #Filter out sequences with stop codons
  if(codon.filter == TRUE){
  x <- codon_filter(x, genetic.code = genetic.code, forward = forward, reverse = reverse)
  }

  #subset to the reading frame
  pos <- get_reading_frame(x, genetic.code = genetic.code, forward = forward, reverse = reverse)

  F_frames <-  as.character(subseq(x, start= pos))

  ent <- vector("list", length=length(F_frames))
  for (l in 1:length(F_frames)){
    ent[[l]] <- c(
      entropy::entropy(table(purrr::map_chr(seq(1, nchar(F_frames[l]), 3), function(i) substr(F_frames[l], i, i))), method=method),
      entropy::entropy(table(purrr::map_chr(seq(2, nchar(F_frames[l]), 3), function(i) substr(F_frames[l], i, i))), method=method),
      entropy::entropy(table(purrr::map_chr(seq(3, nchar(F_frames[l]), 3), function(i) substr(F_frames[l], i, i))), method=method)
      )
    names(ent[[l]]) <- c("pos1", "pos2", "pos3")
  }
  names(ent) <- names(x)
 return(ent)
}
