# Get_ott_taxonomy --------------------------------------------------------

#' Download open tree of life taxonomy
#'
#' @param url a URL to download from, if left blank the latest version will be downloaded
#' @param dest_dir A directory to save the zipped taxonomy database to, if left blank the working directory is selected
#' @param force Whether existing files should be overwritten
#'
#' @return
#' @export
#' @import stringr
#' @importFrom httr GET
#' @importFrom httr write_disk
#' @importFrom utils untar
#' @importFrom rvest html_nodes
#' @importFrom rvest html_attr
#' @importFrom xml2 read_html
#'
download_ott_taxonomy <- function(url, dest_dir, force=FALSE) {

  if (missing(url)) {
    # find the latest version of taxonomy
    download_page <- xml2::read_html("https://tree.opentreeoflife.org/about/taxonomy-version/")
    link_hrefs <- download_page %>%
      rvest::html_nodes("a") %>%
      rvest::html_attr("href")
    url <- grep("http://files.opentreeoflife.org/ott/.*tgz$", link_hrefs, perl = TRUE) %>%
      link_hrefs[.] %>% .[1]
  }

  version <- basename(url) %>%
    stringr::str_remove(".tgz")

  if(missing(dest_dir)){
    message("dest_dir is missing, downloading into /", version)
    dest_dir <- version
  }

  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir) # Create first directory
  }

  # Check if files exist already in dest_dir
  expected_files <- c(
    paste0(dest_dir,"/", version,".tgz" ),
    paste0(dest_dir,"/taxonomy.tsv" ),
    paste0(dest_dir,"/synonyms.tsv" )
  )

  if (any(file.exists(expected_files)) & force == FALSE) {
    message(paste0("Skipped as files already exist in dest_dir, to overwrite set force to TRUE"))
    return(NULL)
  } else  if (any(file.exists(expected_files)) && force == TRUE) {
    unlink(expected_files, recursive = TRUE) # Remove existing version
  }

  dest_file <- file.path(dest_dir, basename(url))
  httr::GET(url, httr::write_disk(dest_file, overwrite=TRUE))

  #unzip file
  utils::untar(dest_file, exdir = ".")
  #Remove download
  file.remove(dest_file)
  message(paste0("Downloaded taxonomy to: ",stringr::str_remove(dest_file,".tgz" ), " \n"))
  return(stringr::str_remove(dest_file,".tgz" ))
}


# get_ott_taxonomy ------------------------------------------------------

#' get_ott_taxonomy
#'
#' @param dir a directory containing the open tree of life taxonomy files obtained from the  `download_ott_taxonomy` function
#' @param quiet Whether progress should be printed to console
#' @param filter_unplaced Whether to filter 'bad' entries. These include
#' incertae_sedis
#' major_rank_conflict
#' unplaced
#' environmental
#' inconsistent
#' extinct
#' hidden
#' hybrid
#' not_otu
#' viral
#' barren
#' See: https://github.com/OpenTreeOfLife/reference-taxonomy/blob/master/doc/taxon-flags.md for more info
#'
#' @return
#' @export
#' @import dplyr
#' @importFrom vroom vroom
#' @importFrom data.table data.table
#'
get_ott_taxonomy <- function(dir=NULL, quiet=FALSE, filter_unplaced=TRUE) {
  if (is.null(dir)){
    input <- NA
    while(!isTRUE(input == "1") && !isTRUE(input == "2")) {
      input <- readline(prompt="No directory provided, type '1' if you want to download the ott taxonomy \n")
      if(input == "1") {
        dir <- download_ott_taxonomy()
      } else {stop("Stopped function")}
    }
  }
  if(!quiet){message("Building data frame\n")}
  file <- normalizePath(paste0(dir, "/taxonomy.tsv"))

  if (filter_unplaced==TRUE){
    remap <- vroom::vroom(file, delim="\t|\t") %>%
      dplyr::filter(!grepl("incertae_sedis|incertae_sedis$|major_rank_conflict|unplaced|environmental|inconsistent|extinct|hidden|hybrid|not_otu|viral|barren", flags)) %>%
      dplyr::mutate(rank = stringr::str_replace(rank, pattern="no rank - terminal", replacement="terminal")) %>%
      dplyr::select(uid, parent_uid, name, rank, sourceinfo, flags) %>%
      dplyr::rename(tax_id = uid, parent_taxid = parent_uid, tax_name = name)
  }else {
    remap <- vroom::vroom(file, delim="\t|\t") %>%
      dplyr::select(uid, parent_uid, name, rank, sourceinfo, flags) %>%
      dplyr::rename(tax_id = uid, parent_taxid = parent_uid, tax_name = name)
  }
  d.dt <- data.table(remap, key="tax_id")
  db <- d.dt[, list(sourceinfo = unlist(strsplit(sourceinfo, ",")), tax_name, parent_taxid, rank, flags), by=tax_id
             ][, c("source", "id") := tstrsplit(sourceinfo, ":", fixed=TRUE)
               ][,c('sourceinfo') :=  .(NULL)]
  attr(db,'type') <- 'OTT'
  return(db)
}

# map_to_ott  ------------------------------------------------------------

#' Map taxa to open tree of life
#'
#' @param x a DNAbin
#' @param db an OTT taxonomic database
#' @param from the existing taxonomic ID format
#' @param resolve_synonyms Whether to resolve synonyms
#' @param dir A directory containing the OTT taxonomy, required if resolve synonyms is true
#' @param filter_unplaced Whether to filter 'bad' entries. These include
#' incertae_sedis
#' major_rank_conflict
#' unplaced
#' environmental
#' inconsistent
#' extinct
#' hidden
#' hybrid
#' not_otu
#' viral
#' barren
#' See: https://github.com/OpenTreeOfLife/reference-taxonomy/blob/master/doc/taxon-flags.md for more info
#' @param remove_na Whether taxa that could not be mapped to the open tree of life should be removed
#' @param quiet Whether progress should be printed to console
#'
#' @return
#' @export
#' @import dplyr
#' @import stringr
#' @importFrom tidyr separate
#' @importFrom tibble as_tibble
#' @importFrom ape as.DNAbin
#' @importFrom methods is
#'
map_to_ott <- function(x, db, from="ncbi", resolve_synonyms=TRUE, dir=NULL, filter_unplaced=TRUE, remove_na = FALSE, quiet=FALSE){
  time <- Sys.time() # get time

  #input checks
  if(resolve_synonyms && is.null(dir)){
    stop("If resolve_synonmys is true, a directory containing the OTT taxonomy must be provided")
  }
  if(!attr(db,'type')=="OTT") {stop("Error: requires OTT db, generate one with get_ott_taxonomy")}
  if(missing(db) & is.null(dir)) {stop("Error: requires OTT db, generate one with get_ott_taxonomy")}
  if(missing(db) & file.exists(paste0(dir,"/taxonomy.tsv"))){
    message("No db provided but dir provided, creating new db")
    db <- get_ott_taxonomy(dir=dir, filter_unplaced = filter_unplaced)
  }
  if (!from %in% unique(db$source)){ stop("Error: 'from' is not in db")}

  #Check input format
  if (methods::is(x, "DNAbin")) {
    message("Input is DNAbin")
    tax <- names(x) %>%
      stringr::str_split_fixed(";", n = 2) %>%
      tibble::as_tibble() %>%
      tidyr::separate(col = V1, into = c("acc", "id"), sep = "\\|") %>%
      dplyr::rename(tax_name = V2)
  } else if (methods::is(x, "DNAStringSet")) {
    message("Input is DNAStringSet")
    tax <- ape::as.DNAbin(x) %>%
      names() %>%
      stringr::str_split_fixed(";", n = 2) %>%
      tibble::as_tibble() %>%
      tidyr::separate(col = V1, into = c("acc", "id"), sep = "\\|") %>%
      dplyr::rename(tax_name = V2)
  }else  if (methods::is(x, "character") && (str_detect(x, "\\|") & str_detect(x, ";"))) {
    message("Detected | and ; delimiters, assuming 'Accession|taxid;Genus Species' format")
    tax <- x %>%
      stringr::str_split_fixed(";", n = 2) %>%
      tibble::as_tibble() %>%
      tidyr::separate(col = V1, into = c("acc", "id"), sep = "\\|") %>%
      dplyr::rename(tax_name = V2)
  } else  if (methods::is(x, "character") && !(stringr::str_detect(x, "\\|") && stringr::str_detect(x, ";"))) {
    message("Did not detect | and ; delimiters, assuming a vector of species names")
    tax <- data.frame(acc = as.character(NA), id=as.character(NA), tax_name = x, stringsAsFactors = FALSE)
  }else (stop("x must be DNA bin or character vector"))

  #Read in synonyms DB
  #TODO: make synonyms df an attribute of the DB object
  if(resolve_synonyms == TRUE){ syn <- parse_ott_synonyms(dir=dir)}

  #Get tax
  tax <- tax %>%
    dplyr::left_join (db %>%   # First map by id
                        dplyr::filter(source==!!from) %>%
                        dplyr::select(-source) %>%
                        dplyr::rename(tax_name.x = tax_name), by = "id") %>%
    dplyr::left_join (db %>%   # Next map by tax_name
                        dplyr::select(-id, -source ) %>%
                        dplyr::filter(!duplicated(tax_name)), by = "tax_name") # then map by name

  #Resolve synonyms
  if(resolve_synonyms == TRUE){
    tax <- tax %>%
      dplyr::left_join(syn %>%
                         dplyr::rename(tax_name.z = tax_name, tax_name = synonym, tax_id.z = tax_id) %>%
                         dplyr::filter(!duplicated(tax_name)),
                       by = "tax_name") %>%
      dplyr::mutate(tax_id = dplyr::case_when(
        !is.na(tax_id.z) ~ tax_id.z, #If synonym was found, use it
        !is.na(tax_id.x) & is.na(tax_id.z) ~ tax_id.x, #If no synonym was found, but an ID match was found use it
        is.na(tax_id.x) & is.na(tax_id.z) & !is.na(tax_id.y)  ~ tax_id.y #If no synonym and no ID match found, but a name match was, use it
      ),
      tax_name = dplyr::case_when(
        !is.na(tax_name.z) ~ tax_name.z, #If synonym was found, use it
        !is.na(tax_name.x) & is.na(tax_name.z) ~ tax_name.x, #If no synonym was found, but an ID match was found use it
        is.na(tax_name.x) & is.na(tax_name.z) ~ tax_name  #If no synonym was found, and no ID match, retain current name
      ))

    if (filter_unplaced == TRUE){ #ensure resolving synonyms didnt introduce bads
      bads <- db %>%
        dplyr::filter(grepl("incertae_sedis,|incertae_sedis$|major_rank_conflict|unplaced|environmental|inconsistent|extinct|hidden|hybrid|not_otu|viral|barren", flags))
      tax <- tax %>%
        dplyr::mutate(tax_id = dplyr::case_when( #Ensure no
          tax_name %in% bads$name ~  as.numeric(NA),
          !tax_name %in% bads$name ~ tax_id
        )) %>%
        dplyr::mutate(name = paste0(acc,"|", tax_id,";",tax_name))
    } else if (filter_unplaced == FALSE){
      tax <- tax %>%
        dplyr::mutate(name = paste0(acc,"|", tax_id,";",tax_name))
    }
  } else if( resolve_synonyms == FALSE){
    tax <- tax %>%
      dplyr::mutate(tax_id = dplyr::case_when(
        !is.na(tax_id.x) ~ tax_id.x,
        is.na(tax_id.x) & !is.na(tax_id.y) ~ tax_id.y
      ),
      tax_name = dplyr::case_when(
        is.na(tax_name.x) ~ tax_name,
        !is.na(tax_name.x) ~ tax_name.x
      )) %>%
      dplyr::mutate(name = paste0(acc,"|", tax_id,";",tax_name))
  }

  #Replace names
  if (methods::is(x, "DNAbin") | methods::is(x, "DNAStringSet")) {
    names(x) <- tax$name
  } else  if (methods::is(x, "character")) {
    x <- tax$name
  }
  time <- Sys.time() - time
  if (!quiet) (message(paste0("Translated ",  length(x)," tax_ids from ",from, " to Open tree of life in ", format(time, digits = 2))))

  # Filter NA's
  if (remove_na ==TRUE){
    remove <- tax %>%
      dplyr::filter(is.na(tax_id)) %>%
      dplyr::mutate(name = paste0(acc,"|", tax_id,";",tax_name))

    if(methods::is(x, "DNAbin") | methods::is(x, "DNAStringSet")){
      x[names(x) %in% remove$name] <- NULL
    }else if (methods::is(x, "character")){
      x[x %in% remove$name] <- NULL
    }
    if(!quiet){message(paste0("Removed ", nrow(remove), " sequences that could not be mapped to OTT\n"))}
  }
  return(x)
}

# Parse Synonyms  -----------------------------------------------------------

#' parse the open tree of life synonyms file
#'
#' @param dir a directory containing the open tree of life taxonomy files obtained from the  `download_ott_taxonomy` function
#' @param quiet Whether progress should be printed to console
#'
#' @return
#' @export
#' @import stringr
#' @import dplyr
#' @importFrom vroom vroom
#'
#'
parse_ott_synonyms <- function(dir=NULL, quiet=FALSE) {
  if (is.null(dir)){
    stop("ERROR: provide a directory containing ott taxonomy")
  }
  if(!quiet){message("Building synonyms data frame\n")}
  file <- normalizePath(paste0(dir, "/synonyms.tsv"))
  out <- vroom::vroom(file, delim="\t|\t" )%>%
    dplyr::mutate(tax_name = uniqname %>%
                    stringr::str_remove(pattern="^.*for ")%>%
                    stringr::str_remove(pattern="\\).*$")%>%
                    stringr::str_remove(pattern="\\(.*$")
    ) %>%
    dplyr::rename(tax_id = uid, synonym = name) %>%
    dplyr::select(tax_id, tax_name, synonym)
  return(out)
}

# OTT Recursion -----------------------------------------------------------


#' Recursively get lineage from OTT taxid
#' This function derives the full lineage of a taxon ID number from a given taxonomy database
#' @param x A DNAbin, DNAStringSet, taxonomy headers formnatted Acc|taxid;taxonomy, or vector of taxids
#' @param db an OTT taxonomic database
#' @param ranks the taxonomic ranks to filter to. Default is "kingdom", "phylum", "class", "order", "family", "genus", "species"
#' To get strain level ranks, add "terminal" to ranks
#' @param output Options are "tax_name" to return taxonomic names for each rank in the lineage or "tax_id" to return taxonomic ID's.
#' @param multithread Whether multithreading should be used, if TRUE the number of cores will be automatically detected (Maximum available cores - 1), or provide a numeric vector to manually set the number of cores to use. Default is FALSE (single thread)
#' This argument may alternatively be a 'cluster' object, in which case it is the user's responsibility to close the socket connection at the conclusion of the operation,
#' for example by running parallel::stopCluster(cores).
#'
#' @return
#' @export
#' @import parallel
#' @import dplyr
#' @import stringr
#' @importFrom tidyr separate
#' @importFrom tidyr unite
#' @importFrom tidyr pivot_wider
#' @importFrom tibble as_tibble
#' @importFrom methods is
#'
get_ott_lineage <- function(x, db, output="tax_name", ranks = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), multithread = FALSE){
  #Check input format
  if (methods::is(x, "DNAbin")) {
    message("Input is DNAbin")
    lineage <- names(x) %>%
      stringr::str_split_fixed(";", n = 2) %>%
      tibble::as_tibble() %>%
      tidyr::separate(col = V1, into = c("acc", "tax_id"), sep = "\\|") %>%
      dplyr::rename(tax_name = V2)
  } else if (methods::is(x, "DNAStringSet")) {
    message("Input is DNAStringSet")
    lineage <- as.DNAbin(x) %>%
      names() %>%
      stringr::str_split_fixed(";", n = 2) %>%
      tibble::as_tibble() %>%
      tidyr::separate(col = V1, into = c("acc", "tax_id"), sep = "\\|") %>%
      dplyr::rename(tax_name = V2)
  }else  if (methods::is(x, "character") && (stringr::str_detect(x, "\\|") & stringr::str_detect(x, ";"))) {
    message("Detected | and ; delimiters, assuming 'Accession|taxid;Genus Species' format")
    lineage <- x %>%
      stringr::str_split_fixed(";", n = 2) %>%
      tibble::as_tibble() %>%
      tidyr::separate(col = V1, into = c("acc", "tax_id"), sep = "\\|") %>%
      dplyr::rename(tax_name = V2)
  } else  if (methods::is(x, "character") && !(stringr::str_detect(x, "\\|") && stringr::str_detect(x, ";"))) {
    message("Did not detect | and ; delimiters, assuming a vector of tax_id's")
    lineage <- data.frame(acc = as.character(NA), tax_name=as.character(NA), tax_id = x, stringsAsFactors = FALSE)
  }else (stop("x must be DNA bin or character vector"))

  # setup multi threading
  cores <- setup_para(multithread, quiet)

  # Check for duplicated accessions
  if(any(duplicated(lineage$acc))){stop("Duplicated sequence accessions found")}

  #check db
  if(missing(db) | !attr(db,'type')=="OTT") {stop("Error: requires OTT db, generate one with get_ott_taxonomy")}

  tax_ids <- as.numeric(lineage$tax_id)
  db$rank <- as.character(db$rank)
  db$tax_name <- as.character(db$tax_name) # avoid stringsasfactor issues

  #dereplicate to uniques
  uh <- unique(paste(tax_ids))
  pointers <- seq_along(uh)
  names(pointers) <- uh
  pointers <- unname(pointers[paste(tax_ids)])
  tax_ids <- tax_ids[!duplicated(pointers)]

  #Recursive function
  gl1 <- function(tax_id, db, ranks){
    if(is.na(tax_id)) return(NA)
    stopifnot(length(tax_id) == 1 & mode(tax_id) == "numeric")
    res <- character(100)
    resids <- integer(100)
    resnames <- character(100)
    counter <- 1
    index <- match(tax_id, db$tax_id)
    if(is.na(index)){
      # warning(paste("Taxon ID", tax_id, "not found in database\n"))
      return(data.frame(rank= NA,
                        tax_name=NA,
                        tax_id=NA,
                        stringsAsFactors = FALSE))
    }
    repeat{
      if(is.na(index)) break
      # if(length(index) > 1) cat(index, "\n")
      res[counter] <- db$tax_name[index]
      resids[counter] <- db$tax_id[index]
      resnames[counter] <- db$rank[index]
      index <- db$parent_tax_index[index]
      counter <- counter + 1
    }
    #get position of ranks
    pos <- match(ranks, resnames)
    resnames <- resnames[pos]
    res <- res[pos]
    resids <- resids[pos]
    out <- data.frame(rank= ranks,
                      tax_name=res,
                      tax_id=resids,
                      stringsAsFactors = FALSE)
    return(out)
  }
  db$parent_tax_index <- match(db$parent_taxid, db$tax_id)
  ## multithreading
  if(inherits(cores, "cluster")){
    res <- parallel::parLapply(cores, tax_ids, gl1, db, ranks)
  }else if(cores == 1){
    res <- lapply(tax_ids, gl1, db, ranks)
  }else{
    navailcores <- parallel::detectCores()
    if(identical(cores, "autodetect")) cores <- navailcores - 1
    if(!(mode(cores) %in% c("numeric", "integer"))) stop("Invalid 'cores' argument")
    if(cores > 1){
      cl <- parallel::makeCluster(cores)
      res <- parallel::parLapply(cl, tax_ids, gl1, db, ranks)
      parallel::stopCluster(cl)
    }else{
      res <- lapply(tax_ids, gl1, db, ranks)
    }
  }
  res <- res[pointers] #re-replicate
  names(res) <- lineage$acc
  if(output =="tax_name"){
    out <- bind_rows(res, .id="id") %>%
      dplyr::select(-tax_id) %>%
      dplyr::group_by(id) %>%
      tidyr::pivot_wider(
        names_from = rank,
        values_from = tax_name) %>%
      dplyr::ungroup() %>%
      dplyr::bind_cols(lineage) %>%
      tidyr::unite(Acc, c(acc, tax_id), sep = "|") %>%
      dplyr::select(Acc, all_of(ranks), tax_name)
  } else if(output == "tax_id"){
    out <- bind_rows(res, .id="id") %>%
      dplyr::select(-tax_name) %>%
      dplyr::group_by(id) %>%
      tidyr::pivot_wider(
        names_from = rank,
        values_from = tax_id) %>%
      dplyr::ungroup() %>%
      dplyr::bind_cols(lineage) %>%
      tidyr::unite(Acc, c(acc, tax_id), sep = "|") %>%
      dplyr::select(Acc, all_of(ranks), tax_name)
  }
  return(out)
}


# Filter ott ---------------------------------------------------------

#' Filter unplaced taxonomic labels
#' @description
#' Filter flags indicating unplaced taxa in the taxonomic tree. These include:
#' incertae_sedis
#' major_rank_conflict
#' unplaced
#' environmental
#' inconsistent
#' extinct
#' hidden
#' hybrid
#' not_otu
#' viral
#' barren
#' See: https://github.com/OpenTreeOfLife/reference-taxonomy/blob/master/doc/taxon-flags.md for more info
#'
#' @param x A DNAbin, DNAStringSet, taxonomy headers formnatted Acc|taxid;taxonomy, or vector of taxids
#' @param db an OTT taxonomic database
#' @param quiet Whether progress should be printed to the console.
#'
#' @return
#' @export
#' @import stringr
#' @import dplyr
#' @importFrom tibble as_tibble
#' @importFrom tidyr separate
#' @importFrom ape as.DNAbin
#' @importFrom methods is
#'
filter_unplaced <- function(x, db, quiet=FALSE){
  #Check input format
  if (methods::is(x, "DNAbin")) {
    message("Input is DNAbin")
    tax <- names(x) %>%
      stringr::str_split_fixed(";", n = 2) %>%
      tibble::as_tibble() %>%
      tidyr::separate(col = V1, into = c("acc", "tax_id"), sep = "\\|") %>%
      dplyr::rename(tax_name = V2)
  } else if (methods::is(x, "DNAStringSet")) {
    message("Input is DNAStringSet")
    tax <- ape::as.DNAbin(x) %>%
      names() %>%
      stringr::str_split_fixed(";", n = 2) %>%
      tibble::as_tibble() %>%
      tidyr::separate(col = V1, into = c("acc", "tax_id"), sep = "\\|") %>%
      dplyr::rename(tax_name = V2)
  }else  if (methods::is(x, "character") && (str_detect(x, "\\|") & str_detect(x, ";"))) {
    message("Detected | and ; delimiters, assuming 'Accession|taxid;Genus Species' format")
    tax <- x %>%
      stringr::str_split_fixed(";", n = 2) %>%
      tibble::as_tibble() %>%
      tidyr::separate(col = V1, into = c("acc", "tax_id"), sep = "\\|") %>%
      dplyr::rename(tax_name = V2)
  } else  if (methods::is(x, "character") && !(str_detect(x, "\\|") && str_detect(x, ";"))) {
    message("Did not detect | and ; delimiters, assuming a vector of species names")
    tax <- data.frame(acc = as.character(NA), tax_id=as.character(NA), tax_name = x, stringsAsFactors = FALSE)
  }else (stop("x must be DNA bin or character vector"))

  #check db
  if(missing(db) | !attr(db,'type')=="OTT") {stop("Error: requires OTT db, generate one with get_ott_taxonomy")}
  bads <- db %>%
    dplyr::filter(grepl("incertae_sedis,|incertae_sedis$|major_rank_conflict|unplaced|environmental|inconsistent|extinct|hidden|hybrid|not_otu|viral|barren", flags))
  if(nrow(bads)==0){stop("Database is already filtered for bad flags, rerun get_ott_taxonomy() with filter_unplaced=FALSE")}
  tax <- tax %>%
    dplyr::mutate(tax_id = dplyr::case_when( #Ensure no
      tax_id %in% bads$tax_id ~  as.character(NA),
      !tax_id %in% bads$tax_id ~ tax_id
    )) %>%
    dplyr::mutate(name = paste0(acc,"|", tax_id,";",tax_name))

  # Filter NA's
  remove <- tax %>%
    dplyr::filter(is.na(tax_id)) %>%
    dplyr::mutate(name = paste0(acc,"|", tax_id,";",tax_name))

  if(methods::is(x, "DNAbin") | methods::is(x, "DNAStringSet")){
    names(x) <- tax$name
    x <- x[!names(x) %in% remove$name]

  }else if (methods::is(x, "character")){
    x <- tax$name
    x <- x[!x %in% remove$name]
  }
  if(!quiet){message(paste0("Removed ", nrow(remove), " unplaced flagged taxa\n"))}
  return(x)

}

#' Filter infraspecific taxonomic labels
#' @param x A DNAbin, DNAStringSet, taxonomy headers formnatted Acc|taxid;taxonomy, or vector of taxids
#' @param db an OTT taxonomic database
#' @param quiet Whether progress should be printed to the console.
#'
#' @return
#' @export
#' @import stringr
#' @importFrom tidyr separate
#' @import dplyr
#' @importFrom tibble as_tibble
#' @importFrom methods is
#'
filter_infraspecifc <- function(x, db, quiet=FALSE){
  #Check input format
  if (methods::is(x, "DNAbin")) {
    message("Input is DNAbin")
    tax <- names(x) %>%
      stringr::str_split_fixed(";", n = 2) %>%
      tibble::as_tibble() %>%
      tidyr::separate(col = V1, into = c("acc", "tax_id"), sep = "\\|") %>%
      dplyr::rename(tax_name = V2)
  } else if (methods::is(x, "DNAStringSet")) {
    message("Input is DNAStringSet")
    tax <- ape::as.DNAbin(x) %>%
      names() %>%
      stringr::str_split_fixed(";", n = 2) %>%
      tibble::as_tibble() %>%
      tidyr::separate(col = V1, into = c("acc", "tax_id"), sep = "\\|") %>%
      dplyr::rename(tax_name = V2)
  }else  if (methods::is(x, "character") && (str_detect(x, "\\|") & str_detect(x, ";"))) {
    message("Detected | and ; delimiters, assuming 'Accession|taxid;Genus Species' format")
    tax <- x %>%
      stringr::str_split_fixed(";", n = 2) %>%
      tibble::as_tibble() %>%
      tidyr::separate(col = V1, into = c("acc", "tax_id"), sep = "\\|") %>%
      dplyr::rename(tax_name = V2)
  } else  if (methods::is(x, "character") && !(str_detect(x, "\\|") && str_detect(x, ";"))) {
    message("Did not detect | and ; delimiters, assuming a vector of species names")
    tax <- data.frame(acc = as.character(NA), tax_id=as.character(NA), tax_name = x, stringsAsFactors = FALSE)
  }else (stop("x must be DNA bin or character vector"))

  #check db
  if(missing(db) | !attr(db,'type')=="OTT") {stop("Error: requires OTT db, generate one with get_ott_taxonomy")}
  bads <- db %>%
    dplyr::filter(grepl("infraspecific", flags))
  tax <- tax %>%
    dplyr::mutate(tax_id = dplyr::case_when( #Ensure no
      tax_id %in% bads$tax_id ~  as.character(NA),
      !tax_id %in% bads$tax_id ~ tax_id
    )) %>%
    dplyr::mutate(name = paste0(acc,"|", tax_id,";",tax_name))

  # Filter NA's
  remove <- tax %>%
    dplyr::filter(is.na(tax_id)) %>%
    dplyr::mutate(name = paste0(acc,"|", tax_id,";",tax_name))

  if(methods::is(x, "DNAbin") | methods::is(x, "DNAStringSet")){
    names(x) <- tax$name
    x <- x[!names(x) %in% remove$name]

  }else if (methods::is(x, "character")){
    x <- tax$name
    x <- x[!x %in% remove$name]
  }
  if(!quiet){message(paste0("Removed ", nrow(remove), " infraspecific taxa\n"))}

  # check for multiple taxids per name remaining
  checks <- tax %>%
    dplyr::filter(!is.na(tax_id)) %>%
    dplyr::group_by(tax_name) %>%
    dplyr::summarise(count = n_distinct(tax_id))

  if (any(checks$count >1 )) {warning("Multiple tax_ids remain for some tax_names, please check manually")}
  return(x)
}

