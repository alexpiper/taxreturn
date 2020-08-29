#' Download NCBI taxdump
#'
#' @param dest.dir A directory to save the downloaded ncbi taxdump files. If empty a new folder in the working directory will be created
#' @param synonyms Whether synonyms should be included in the database
#' @param force Whether already downloaded files should be overwritten
#'
#' @return
#' @export
#' @import utils
#' @import readr
#'
#' @examples
get_ncbi_lineage <- function(dest.dir, synonyms = TRUE, force=FALSE) {
  if(missing(dest.dir)){
    #message("dest.dir is missing, downloading/looking for files in working directory")
    dest.dir <- getwd()
  }
  tmp <- paste0(dest.dir,"/ncbi_taxdump")
  if (!dir.exists(tmp)) {
    dir.create(tmp) # Create first directory
  }
  if (force == TRUE | !file.exists(paste0(tmp, "/rankedlineage.dmp"))) {
    message("Downloading NCBI taxonomy database")
    fn <- "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
    download.file(fn, destfile = paste0(tmp, "/tmp.tar.gz"))
    message("Extracting data\n")
    test <- utils::untar(tarfile = paste0(tmp, "/tmp.tar.gz"), exdir = tmp)
    if (!identical(test, 0L)) {
      stop(cat(test))
    }
    file.remove(paste0(tmp, "/tmp.tar.gz"))
  }
  message("Building NCBI taxonomy data frame\n")

  ## Read data frame
  #lin <- readr::read_tsv(paste0(tmp, "/rankedlineage.dmp"),
  #                    delim="\t|",
  #                    col_names = c("tax_id", "tax_name", "species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom", "junk"),
  #                    col_types = list(col_integer(),
  #                                     col_character(),
  #                                     col_character(),
  #                                     col_character(),
  #                                     col_character(),
  #                                     col_character(),
  #                                     col_character(),
  #                                     col_character(),
  #                                     col_character(),
  #                                     col_character(),
  #                                     col_character()
#
  #  )
#
  #) %>%
  #  select(-junk)
  lin <- readr::read_tsv(paste0(tmp, "/rankedlineage.dmp"),
                         col_names = c("tax_id", "tax_name", "species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom"),
                         col_types = ("i-c-c-c-c-c-c-c-c-c-")
  )

  # Remove synonyms
  if (synonyms == FALSE) {
    x <- scan(
      file = paste0(tmp, "/names.dmp"), what = "", sep = "\n",
      quiet = TRUE
    )
    syn <- x[grepl("synonym", x)]
    syn <- strsplit(syn, split = "\t")
    syn <- sapply(syn, function(s) s[c(1, 3)])
    syn <- as.data.frame(t(syn), stringsAsFactors = FALSE)
    syn[[1]] <- as.integer(syn[[1]])
    colnames(syn) <- c("taxID", "name")
    lin <- lin %>%
      filter(!tax_id %in% syn[[1]])
  }

  message("Done\n")
  attr(lin,'type') <- 'ncbi'
  return(lin)
}

# NCBI synonyms -----------------------------------------------------------
#' Get NCBI synonyms
#'
#' @param dir A directory containing the NCBI taxonomy that was downloaded using get_ncbi_lineage
#' @param recurse Whether to recurse when searching for dir if dir is NULL.
#' If TRUE recurse fully, if a positive number the number of levels to recurse.
#' @param quiet Whether progress should be printed to console
#'
#' @return
#' @export
#' @import fs
#' @import readr
#' @import dplyr
#'
#' @examples
get_ncbi_synonyms <- function(dir=NULL, recurse=TRUE, quiet=FALSE) {
  if (is.null(dir)){
    if(!quiet){message("Dir is null, finding ncbi_taxdump directory using recursion")}
    dir <- fs::dir_ls(path = getwd(), type = "directory", glob = "*ncbi_taxdump", recurse = recurse)
    if(length(dir) >0){
      if(!quiet){message("ncbi_taxdump found at ", dir)}
    } else{
      stop("ERROR: provide a directory containing ncbi taxonomy")
    }
  }
  if(!quiet){message("Building synonyms data frame\n")}
  file <- normalizePath(paste0(dir, "/names.dmp"))

  #Read in
  splitLines <- do.call(rbind, strsplit(readLines(file),"\t\\|\t?"))
  splitLines<-splitLines[,-3]
  colnames(splitLines) <- c("tax_id","tax_name", "class")
  parsed <- as.data.frame(splitLines)

  out <- parsed %>%
    dplyr::group_by(tax_id) %>%
    dplyr::filter(any(class=="synonym")) %>%
    ungroup() %>%
    dplyr::filter(!class=="scientific name")  %>%
    dplyr::select(-class, synonym = tax_name) %>%
    dplyr::left_join(parsed %>%
                  dplyr::filter(class=="scientific name") %>%
                  dplyr::select(tax_id, tax_name), by="tax_id")

  return(out)
}

# resolve_synonyms_ncbi ---------------------------------------------------

#' resolve_synonyms_ncbi
#'
#' @param x A DNAbin or DNAStringSet Object
#' @param quiet
#'
#' @return
#' @export
#' @import ape
#' @import Biostrings
#' @import stringr
#' @import tibble
#' @import tidyr
#' @import dplyr
#'
#' @examples
resolve_synonyms_ncbi <- function(x, dir=NULL, quiet = FALSE, ...) {
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

  syns <- get_ncbi_synonyms(dir=dir) #...=...

  # if input has sequences, get names
  if (is.seq) {
    query <- names(x) %>%
      stringr::str_split_fixed(";", n = 2) %>%
      tibble::as_tibble() %>%
      tidyr::separate(col = V1, into = c("acc", "tax_id"), sep = "\\|") %>%
      dplyr::rename(query = V2)
  } else if (!is.seq) {
    query <- data.frame(query= x)
    query$tax_id <- "NA" # Add dummy columns
    query$acc <- "NA"
  }

  #Get synonyms to replace
  to_replace <- query %>%
    dplyr::select(-tax_id) %>%
    dplyr::filter(query %in% syns$synonym) %>%
    dplyr::left_join(syns %>% dplyr::rename(query = synonym)) %>%
    dplyr::select(-query) %>%
    dplyr::filter(!duplicated(acc)) #where ar these coming from?

  if(nrow(to_replace) > 0){
  out <- query %>%
    dplyr::rename(tax_name = query) %>%
    dplyr::rows_update(to_replace, by="acc") %>%
    tidyr::unite(col = V1, c("acc", "tax_id"), sep = "|")


  if(is.seq) {
      names(x) <- paste(out$V1, out$tax_name, sep = ";")
    } else if (!is.seq) {
      x  <- paste(out$V1, out$tax_name, sep = ";")
    }
  time <- Sys.time() - time
  if (!quiet) {message(paste0("resolved ", nrow(to_replace), " synonyms in ", format(time, digits = 2)))}
  } else (message("No synonyms detected"))
  return(x)
}

#
#' (Deprecated) Get ranked Lineage
#'
#' @param db
#' @param synonyms
#' @param force
#'
#' @return
#' @export
#'
#' @examples
get_ranked_lineage <- function(db = "NCBI", synonyms = TRUE, force=FALSE) {
  .Deprecated(new="get_ncbi_lineage", old="get_ranked_lineage")
  get_ncbi_lineage(db = "NCBI", synonyms = TRUE, force = FALSE)
}


# ncbi_taxid --------------------------------------------------------------

#' Get ncbi taxid's for a taxon name
#'
#' @param x
#' @param db
#'
#' @return
#' @export
#' @import magrittr
#' @import dplyr
#'
#' @examples
ncbi_taxid <- function(x, db=NULL) {

  if (is.null(db)) { db <- get_ncbi_lineage()}
  out <-  as.data.frame(x) %>%
    magrittr::set_colnames("tax_name") %>%
    dplyr::left_join (db, by = "tax_name") %>%
    dplyr::pull(tax_id)

  return(out)
}


# gid_to_acc --------------------------------------------------------------


#' Convert NCBI gene ids to accession numbers
#'
#' @param ids A character or numeric vector of NCBI gids
#' @param database The origin database of the gids, default 'nuccore'
#' @param chunksize The size of the chunked searches to conduct.
#' Warning, chunk sizes over 300 can be too big for the NCBI servers.
#' @param multithread Whether multithreading should be used
#' @param quiet (Optional) Print text outputdatabase
#'
#' @return
#' @export
#' @import future
#' @import furrr
#' @import rentrez
#' @import purrr
#'
#'
#' @examples
gid_to_acc <- function(ids, database="nuccore", chunksize=300, multithread=TRUE, progress=FALSE, quiet=FALSE){
  if(!class(ids) %in% c("character", "numeric")){
    stop("input ids must be a character or numeric vector of NCBI gids")
  }
  time <- Sys.time() # get time

  # split search into chunks
  chunks <- split(ids, ceiling(seq_along(ids)/chunksize))
  if(!quiet){message(paste0("Converting ", length(ids), " NCBI gids to accession numbers in ", length(chunks), " chunks"))}

  # setup multithreading
  ncores <- future::availableCores() -1
  if(isTRUE(multithread)){
    cores <- ncores
    if(!quiet){message("Multithreading with ", cores, " cores")}
    future::plan(future::multiprocess, workers=cores)
  } else if (is.numeric(multithread) & multithread > 1){
    cores <- multithread
    if(cores > ncores){
      cores <- ncores
      message("Warning: the value provided to multithread is higher than the number of cores, using ", cores, " cores instead")
    }
    if(!quiet){message("Multithreading with ", cores, " cores")}
    future::plan(future::multiprocess, workers=cores)
  } else if(isFALSE(multithread) | multithread==1){
    future::plan(future::sequential)
  } else (
    stop("Multithread must be a logical or numeric vector of the numbers of cores to use")
  )

  #Main function
  out <- furrr::future_map(chunks, function(x){
    upload <- rentrez::entrez_post(db=database, id=x)
    dl <- rentrez::entrez_fetch(db = database, web_history = upload, rettype = "acc", retmax = chunksize)
    acc <- readLines(textConnection(dl))
    acc <- acc[!acc==""]
    return(acc)
  }, .progress = progress) %>%
    unlist(use.names=FALSE)

  #finished
  time <- Sys.time() - time
  if (!quiet) {message(paste0("Sucessfully converted ", length(out), " of ", length(ids), " NCBI gids to accession numbers in ", format(time, digits = 2)))}
  return(out)
}
