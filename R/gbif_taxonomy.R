### Get GBIF Taxonomy

# modified from traitdataform pacakge https://github.com/EcologicalTraitData/traitdataform/issues/35

#' Resolve GBIF
#'
#' @param x
#' @param subspecies
#' @param higherrank
#' @param verbose
#' @param fuzzy
#' @param conf_threshold
#' @param resolve_taxonomy
#'
#' @import taxize
#' @import purrr
#' @return
#'
#' @examples
resolve_gbif <- function(x, subspecies = TRUE, higherrank = TRUE, verbose = FALSE,
                         fuzzy = TRUE, conf_threshold = 90, resolve_taxonomy = TRUE) {
  matchtype <- status <- confidence <- NULL

  ## Get GBIF data - this needs to be sped up majorly, would be nice to move to taxize::db
  # temp <- taxize::get_gbifid_(x, messages = verbose)

  temp <- x %>%
    purrr::map(safely(taxize::get_gbifid_)) %>%
    purrr::map("result") %>%
    purrr::flatten()

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
      temp[[i]] <- subset(temp[[i]], temp[[i]]$confidence == max(temp[[i]]$confidence))
      if (nrow(temp[[i]]) > 1) {
        temp[[i]] <- temp[[i]][1, ]
        warning_i <- paste(warning_i, "Selected first of multiple equally ranked concepts!")
      }
    }
    if (!any(temp[[i]]$status == "ACCEPTED") & any(temp[[i]]$status == "SYNONYM")) {
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
          temp[[i]] <- subset(temp[[i]], temp[[i]]$confidence ==  max(temp[[i]]$confidence))
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
                           sep = ""  )
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
        warning_i <- paste(warning_i, "No matching species concept! Entry has been mapped to higher taxonomic rank.")
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
        author = sub(paste0( temp[[i]]$canonicalname,   " "  ), "", temp[[i]]$scientificname), taxonRank = temp[[i]]$rank,
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

# Resolve taxonomy -------------------------------------------------------

#' Resolve taxonomic synonyms
#' @description This function takes a DNAbin object, or a list of species and uses the Global Biodiversity Information Facility (GBIF) to resolve taxonomic synonyms
#'
#' @param x A DNAbin or DNAStringset object
#' @param subspecies Whether subspecies should be included
#' @param quiet Whether progress should be printed to the console.
#' @param missing How to handle cases where the updated taxonomic name does not have a taxonomic ID in the NCBI database.
#' Options include "ignore", which will update name but keep old taxid, "keepold" which will keep the old taxname and taxid,
#' and "remove" which will remove all synonyms that dont have a taxid in the NCBI database.
#' @param higherrank Whether taxa that were not found should be mapped to a higher rank
#' @param fuzzy Whether fuzzy matching of taxa names should be used
#'
#' @return This returns a DNAbin with renamed taxa
#' @export
#'
#' @import ape
#' @import Biostrings
#' @import stringr
#' @import tibble
#' @import tidyr
#' @import dplyr
#' @import insect
#'
#' @examples
resolve_taxonomy <- function(x, subspecies = FALSE, quiet = TRUE, missing = "ignore", higherrank = FALSE, fuzzy = TRUE) {
  time <- Sys.time() # get time

  if (quiet == TRUE) {
    verbose <- FALSE
  } else {
    (verbose <- TRUE)
  }
  db <- get_ncbi_lineage(synonyms = TRUE, force=FALSE)

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
    dplyr::left_join(db %>% select(tax_name, tax_id) %>% rename(scientificNameStd = tax_name), by="scientificNameStd") %>%
    dplyr::rename(taxidnew = tax_id) %>%
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
        )) %>%
        dplyr::mutate(taxid = case_when(
          !is.na(taxidnew) ~ as.character(taxidnew),
          is.na(taxidnew) ~ taxid
        )) %>%
        dplyr::select(acc, taxid, query)
    } else if (missing == "keepold") { # In cases where the updated synonym does not have a taxonomic ID in the NCBI database, keep old taxname

      query <- query %>%
        dplyr::left_join(out, by = "query") %>%
        dplyr::mutate(query = case_when(
          !is.na(scientificNameStd) & !is.na(taxidnew) ~ scientificNameStd,
          TRUE ~ query
        )) %>% # Catch all for anything that is not above case
        dplyr::mutate(taxid = case_when(
          !is.na(scientificNameStd) & !is.na(taxidnew) ~ as.character(taxidnew),
          is.na(taxidnew) ~ taxid
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
        #dplyr::filter(keepcol == TRUE) %>%
        dplyr::mutate(query = case_when(
          keepcol ==FALSE ~ "REMOVE",
          !is.na(scientificNameStd) & keepcol ==TRUE ~ scientificNameStd,
          TRUE ~ query
        )) %>% # Catch all for anything that is not above case
        dplyr::mutate(taxid = case_when(
          !is.na(taxidnew) ~ as.character(taxidnew),
          is.na(taxidnew) ~ taxid
        ))
      removed <- sum(query$keepcol==FALSE)

      query <- query %>%
        dplyr::select(acc, taxid, query)
    }

    if(is.seq == TRUE) {
      query <- query %>%
        tidyr::unite(col = V1, c("acc", "taxid"), sep = "|")

      names(x) <- paste(query$V1, query$query, sep = ";")
      x <- insect::subset.DNAbin(x, subset = !str_detect(names(x), "REMOVE"))
    } else if (is.seq == FALSE) {
      x <- query$query[which(!str_detect(query$query, "REMOVE"))]

    }

  } else (message("No synonyms detected"))

  time <- Sys.time() - time
  if (!quiet & !missing == "remove") {message(paste0("resolved ", nrow(out), " synonyms in ", format(time, digits = 2)))}
  if (!quiet & missing == "remove") {message(paste0("resolved ", nrow(out), " synonyms in ", format(time, digits = 2)," and removed ", removed, " taxa with no updated synonym in NCBI removed as missing = 'remove'"))}
  return(x)
}
