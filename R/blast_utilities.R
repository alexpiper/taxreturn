
# .findExecutable ---------------------------------------------------------


#' Find executable
#'
#' @param exe the name of the executable to find
#' @param quiet Whether errors should be printed to console
#'
#' @return
#' @export
#'
.findExecutable <- function(exe, quiet=FALSE) {
  path <- Sys.which(exe)
  if(all(path=="")) {
    if(!quiet) stop("Executable for ", paste(exe, collapse=" or "), " not found! Please make sure that the software is correctly installed and, if necessary, path variables are set.", call.=FALSE)
    return(character(0))
  }

  path[which(path!="")[1]]
}


# Install BLAST -----------------------------------------------------------

#' Install BLAST
#'
#' @param url (Optional) Default will search for the latest version
#' URL to retrieve BLAST version from.
#' @param dest_dir (Optional)  Default "bin"
#' Directory to install BLAST within.
#' @param force Whether existing installs should be forcefully overwritten
#'
#' @return
#' @export
#' @import stringr
#' @importFrom RCurl getURL
#' @importFrom httr GET
#' @importFrom httr write_disk
#' @importFrom utils untar
#'
blast_install <- function(url, dest_dir = "bin", force = FALSE) {

  # get start time
  time <- Sys.time()
  # get OS
  localos <- Sys.info()["sysname"]

  if (missing(url)) {
    # find the latest version of BLAST
    url <- 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/'
    filenames <- RCurl::getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE) %>%
      stringr::str_split("\r*\n") %>%
      unlist()

    if(localos == "Windows"){
    url <- filenames[stringr::str_detect(filenames,"win64.tar.gz$")] %>%
      paste0(url,.)
    } else if(localos == "Darwin"){
      url <- filenames[stringr::str_detect(filenames,"macosx.tar.gz$")] %>%
        paste0(url,.)
    } else if(localos == "unix"){
      url <- filenames[stringr::str_detect(filenames,"linux.tar.gz$")] %>%
        paste0(url,.)
    }

  }

  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir) # Create first directory
  }

  blast_version <- basename(url) %>% stringr::str_replace("(-x64)(.*?)(?=$)", "")
  if (dir.exists(paste0(dest_dir, "/", blast_version)) && force == FALSE) {
    message("Skipped as BLAST already exists in directory, to overwrite set force to TRUE")
    return(NULL)
  } else  if (dir.exists(paste0(dest_dir, "/", blast_version)) && force == TRUE) {
    unlink(paste0(dest_dir, "/", blast_version), recursive = TRUE) # Remove old version
  }

  destfile <- file.path(dest_dir, basename(url))
  if (exists(destfile)) {
    file.remove(destfile) # Remove old zip file
  }

  #Download and unzip
  httr::GET(url, httr::write_disk(destfile, overwrite=TRUE))
  utils::untar(destfile, exdir = dest_dir)
  file.remove(destfile)

  #Set new $Paths variable for mac & linux
  if(localos == "Darwin" | localos == "unix"){
    old_path <- Sys.getenv("PATH")
    install_path <- list.dirs(dest_dir, full.names = TRUE)[str_detect(list.dirs(dest_dir, full.names = TRUE),"/bin$")]
    Sys.setenv(PATH = paste(old_path, normalizePath(install_path), sep = ":"))
  }

  time <- Sys.time() - time
  message(paste0("Downloaded ", blast_version, " in ", format(time, digits = 2)))
}


# Make Blast DB -----------------------------------------------------------

#' Make blast Database
#'
#' @param file (Required) A fasta file to create a database from.
#' @param dbtype (Optional) Molecule type of database, accepts "nucl" for nucleotide or "prot" for protein.
#' @param args (Optional) Extra arguments passed to BLAST
#' @param quiet (Optional) Whether progress should be printed to console, default is FALSE
#' @param remove_gaps (Optional) Whether gaps should be removed from the fasta file. Note that makeblastdb can fail if there are too many gaps in the sequence.
#'
#' @return
#' @export
#' @import stringr
#' @import ape
#' @importFrom R.utils gunzip
#'
make_blast_db <- function (file, dbtype = "nucl", args = NULL, quiet = FALSE, remove_gaps=TRUE) {
  time <- Sys.time() # get time
  .findExecutable("makeblastdb") # Check blast is installed
  if (is.null(args)){args <- ""}

  # Create a temp file for new blast DB
  tmpfile <- tempfile()

  if(remove_gaps) {
    seqs <- ape::del.gaps(ape::read.FASTA(file))
    write_fasta(seqs, tmpfile, compress = FALSE)
  } else {
    if (stringr::str_detect(file, ".gz")) {
      message("Unzipping file")
      R.utils::gunzip(filename=file, destname=tmpfile, remove=FALSE, overwrite=TRUE)
    } else {
      file.copy(file, tmpfile)
    }
  }

  # Run makeblastdb executable
  results <- system2(command = .findExecutable("makeblastdb"),
                     args = c("-in", tmpfile, "-dbtype", dbtype, args),
                     wait = TRUE,
                     stdout = TRUE)
  time <- Sys.time() - time
  if (!quiet) (message(paste0("made BLAST DB in ", format(time, digits = 2))))
  return(tmpfile)
}


#' Show BLAST parameters
#'
#' @param type (Required) Which BLAST function to display help page for
#'
#' @return
#' @export
#'
blast_params <- function(type = "blastn") {
  system(paste(.findExecutable(c(type)), "-help"))
}


# BLAST -------------------------------------------------------------------

#' Run BLAST search
#'
#' @param query (Required) Query sequence. Accepts a DNABin object, DNAStringSet object, Character string, or filepath.
#' @param db (Required) Reference sequences to conduct search against. Accepts a DNABin object, DNAStringSet object, Character string, or filepath.
#' If DNAbin, DNAStringSet or character string is provided, a temporary fasta file is used to construct BLAST database.
#' If db is set to "remote", this will conduct a search against NCBI nucleotide database.
#' @param type (Required) type of search to conduct, default 'blastn'
#' @param evalue (Required) Minimum evalue from search
#' @param output_format The output format to be returned.
#'  Default is tabular, which returns 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcov
#' See https://www.ncbi.nlm.nih.gov/books/NBK279684/ for custom formatting information
#' @param args (Optional) Extra arguments passed to BLAST
#' @param ungapped Whether ungapped alignment should be conducted. Default is FALSE.
#' @param quiet (Optional) Whether progress should be printed to console, default is FALSE
#' @param multithread Whether multithreading should be used, if TRUE the number of cores will be automatically detected, or provided a numeric vector to manually set the number of cores to use
#' @param remove_db_gaps Whether gaps should be removed from the fasta file used for the database. Note that makeblastdb can fail if there are too many gaps in the sequence.
#'
#' @return
#' @export
#' @import stringr
#' @import dplyr
#' @import future
#' @importFrom tibble enframe
#' @importFrom Biostrings writeXStringSet
#' @importFrom Biostrings DNA_ALPHABET
#' @importFrom ape as.DNAbin
#' @importFrom ape del.gaps
#' @importFrom stringr str_remove
#' @importFrom stringr str_extract
#' @importFrom stringr str_split_fixed
#' @importFrom stringr str_to_upper
#' @importFrom stringr str_split
#' @importFrom future availableCores
blast <- function (query, db, type="blastn", evalue = 1e-6,
                   output_format = "tabular", args=NULL, ungapped=FALSE,
                   quiet=FALSE, multithread=FALSE, remove_db_gaps = TRUE){
  time <- Sys.time() # get time
  # Create temp files
  tmp <- tempdir()
  tmpquery <- paste0(tmp, "/tmpquery.fa")
  tmpdb <- paste0(tmp, "/tmpquery.fa")

  .findExecutable(type) # check blast is installed

  #Setup multithreading
  ncores <- future::availableCores() -1
  if((isTRUE(multithread) | is.numeric(multithread) & multithread > 1) & db=="remote"){
    stop("Multithreading must be set to false for remote searches")
  } else if(isTRUE(multithread)){
    cores <- ncores
    if(!quiet){message("Multithreading with ", cores, " cores")}
  } else if (is.numeric(multithread) & multithread > 1){
    cores <- multithread
    if(cores > ncores){
      cores <- ncores
      message("Warning: the value provided to multithread is higher than the number of cores, using ", cores, " cores instead")
    }
    if(!quiet){message("Multithreading with ", cores, " cores")}
  } else if(isFALSE(multithread) | multithread==1){
    cores <- 1
  } else (
    stop("Multithread must be a logical or numeric vector of the numbers of cores to use")
  )
  nthreads <- ifelse(cores > 1, paste0("-num_threads ", unname(cores)), "")

  # Check outfmt
  if(output_format=="tabular"){
    #Custom tabular format
    parsecols <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                   "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs", "qcovhsp")
    outfmt <- paste0("\"",6," ", paste(parsecols, collapse=" "),"\"")
  } else if(is.numeric(output_format)){
    outfmt <- output_format
    parsecols <- NULL
  } else if (!is.na(stringr::str_extract(output_format, "."))){
    outfmt <- paste0("\"",output_format,"\"")
    parsecols <- output_format %>%
      stringr::str_remove("^..") %>%
      stringr::str_split_fixed("\ ", n=Inf) %>%
      as.character()
  }

  # Database
  if(db=="remote"){
    db <- "nt"
    remote <- "-remote"
  } else if(inherits(db, "DNAbin")){
    if (!quiet) { message("Database input is DNAbin: Creating temporary blast database") }
    write_fasta(db, tmpdb)
    db <- make_blast_db(tmpdb, remove_gaps = remove_db_gaps)
    remote <- ""
  } else if (inherits(db, "DNAString") | inherits(db, "DNAStringSet")){
    if (!quiet) { message("Database input is DNAStringSet: Creating temporary blast database") }
    Biostrings::writeXStringSet(db, tmpdb)
    db <- make_blast_db(tmpdb, remove_gaps = remove_db_gaps)
    remote <- ""
  } else if (inherits(db, "character") &&  all(stringr::str_to_upper(stringr::str_split(db,"")[[1]]) %in% Biostrings::DNA_ALPHABET)) { # Handle text input
    if (!quiet) { message("Database input is character string: Creating temporary blast database") }
    if (nchar(db[1]) == 1) {db <- paste0(db, collapse = "")}
    db <- char2DNAbin(db)
    write_fasta(db, tmpdb)
    db <- make_blast_db(tmpdb, remove_gaps = remove_db_gaps)
    remote <- ""
  } else if (inherits(db, "character") &&  file.exists(file.path(db))){ # Handle filename
    db <- make_blast_db(db, remove_gaps = remove_db_gaps)
    remote <- ""
  } else {
    stop("Invalid BLAST database")
  }

  # Query
  if(inherits(query, "DNAbin")){
    if (!quiet) { message("Query is DNAbin: Creating temporary fasta file") }
    query <- ape::del.gaps(query)
    write_fasta(query, tmpquery)
    input <- tmpquery
  } else if (inherits(query, "DNAString") | inherits(query, "DNAStringSet")){
    if (!quiet) { message("Query is DNAString: Creating temporary fasta file") }
    query <- ape::as.DNAbin(query)
    query <- ape::del.gaps(query)
    write_fasta(query, tmpquery)
    input <- tmpquery
  }else if (inherits(query, "character") &&  all(stringr::str_to_upper(stringr::str_split(query,"")[[1]]) %in% Biostrings::DNA_ALPHABET)) { # Handle text query
    if (!quiet) { message("Query is character string: Creating temporary fasta file") }
    if (nchar(query[1]) == 1) {query <- paste0(query, collapse = "")}
    query <- char2DNAbin(query)
    query <- ape::del.gaps(query)
    write_fasta(query, tmpquery)
    input <- tmpquery
  } else if (inherits(query, "character") &&  file.exists(file.path(query))){ # Handle filenames
    input <- query
  }else {
    stop("Invalid BLAST query")
  }

  #Setup ungapped
  ungapped <- ifelse(ungapped, "-ungapped ","")

  #  Conduct BLAST search
  if (!quiet) { message("Running BLASTn query: ", paste(c("-db", db,
                                                          "-query", input,
                                                          remote,
                                                          "-outfmt", outfmt,
                                                          "-evalue", evalue,
                                                          args,
                                                          nthreads,
                                                          ungapped), collapse=" "))}
  results <- system2(command = .findExecutable(type),
                     args = c("-db", db,
                              "-query", input,
                              remote,
                              "-outfmt", outfmt,
                              "-evalue", evalue,
                              args,
                              nthreads,
                              ungapped),
                     wait = TRUE,
                     stdout = TRUE)

  # Check for errors and stop if detected
  if(any(stringr::str_detect(results, "Error:"))){
    err_to_print <- results[stringr::str_detect(results, "Error:")]
    stop(err_to_print[1])
  }

  # Remove the error message about nucleotides in title
  results <- results[!stringr::str_detect(results, "Title ends with at least 20 valid nucleotide characters")]

  # Parse BLAST results
  if(!is.null(parsecols)){
    out <- results %>%
      tibble::enframe() %>%
      tidyr::separate(col = value, into = parsecols,  sep = "\t", convert = TRUE)
  } else{
    message("Warning, result parsing is currently only supported for output_format = 'tabular', returning raw results")
    out <- results %>%
      tibble::enframe()
  }
  time <- Sys.time() - time
  if (!quiet) (message(paste0("finished BLAST in ", format(time, digits = 2))))

  # Clean up files
  if(file.exists(tmpdb)){file.remove(tmpdb)}
  if(file.exists(tmpquery)){file.remove(tmpquery)}
  return(out)
}


# BLAST_top_hit -----------------------------------------------------------

#' BLAST Top Hit
#'
#' @description Conduct BLAST search and return top hit
#' @param query (Required) Query sequence. Accepts a DNABin object, DNAStringSet object, Character string, or filepath.
#' @param db (Required) Reference sequences to conduct search against. Accepts a DNABin object, DNAStringSet object, Character string, or filepath.
#' If DNAbin, DNAStringSet or character string is provided, a temporary fasta file is used to construct BLAST database
#' @param type (Required) type of search to conduct, default 'blastn'
#' @param identity (Required) Minimum percent identity cutoff. Note that this is calculated using all alignments for each query-subject match.
#' @param coverage (Required) Minimum percent query coverage cutoff. Note that this is calculated using all alignments for each query-subject match.
#' @param evalue (Required) Minimum expect value (E) for saving hits
#' @param max_target_seqs (Required) Number of aligned sequences to keep. Even if you are only looking for 1 top hit keep this higher for calculations to perform properly.
#' @param max_hsp (Required) Maximum number of HSPs (alignments) to keep for any single query-subject pair.
#' @param ranks (Required) The taxonomic ranks contained in the fasta headers
#' @param delim (Required) The delimiter between taxonomic ranks in fasta headers
#' @param tie How to handle ties in top hit results. Options are to break ties by selecting the first hit (Default), or return all tied hits.
#' @param args (Optional) Extra arguments passed to BLAST
#' @param quiet (Optional) Whether progress should be printed to console, default is FALSE
#' @param remove_db_gaps Whether gaps should be removed from the fasta file used for the database. Note that makeblastdb can fail if there are too many gaps in the sequence.
#'
#' @return
#' @export
#' @import dplyr
#' @importFrom tidyr separate
#'
blast_top_hit <- function(query, db, type="blastn",
                          identity=95, coverage=95, evalue=1e06, max_target_seqs=5, max_hsp=5,
                          ranks=c("Kingdom", "Phylum","Class", "Order", "Family", "Genus", "Species"), delim=";",
                          resolve_ties="first", args=NULL, remove_db_gaps = TRUE, quiet=FALSE,...){

  # set input filters in advance to speed up blast
  args <- paste("-perc_identity", identity, "-max_target_seqs", max_target_seqs, "-max_hsps", max_hsp, args)


  #Conduct BLAST
  result <- blast(query=query, type=type, db=db,
                  evalue = evalue,
                  args=args,
                  output_format = '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend qlen sstart send evalue bitscore qcovs',
                  remove_db_gaps = remove_db_gaps) %>%
    dplyr::filter(!is.na(sseqid))

  # Check if ranks set up correctly
  if(!length(unlist(str_split(result$stitle[1], ";"))) == length(c("acc", ranks))){
    stop("number of ranks does not match database!")
  }

  #Subset to top hit
  # Full-length pid from: https://github.com/McMahonLab/TaxAss/blob/master/tax-scripts/calc_full_length_pident.R
  top_hit <- result %>%
    dplyr::mutate( q_align = qend - qstart + 1) %>% # Calculate number of aligned bases per hsp -  add 1 b/c start is 1 instead of 0.
    dplyr::mutate(full_pident = (pident * length) / (length - q_align + qlen) ) %>% # Calculate full pid for each hsp
    dplyr::group_by(qseqid, sseqid, stitle) %>%
    dplyr::summarise(pident = sum(full_pident),# Summarise to per subject-query hits
              qcovs = unique(qcovs), #qcovs are calculated per subjec talready
              max_score = max(bitscore), # Calculated as per NCBI web blast
              total_score = sum(bitscore), # Calculated as per NCBI web blast
              evalue = min(evalue) # Calculated as per NCBI web blast
    )  %>%
    dplyr::filter(pident > identity, qcovs > coverage) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(qseqid) %>%
    dplyr::top_n(1,total_score) %>%
    dplyr::top_n(1,max_score) %>%
    dplyr::top_n(1,qcovs) %>%
    dplyr::top_n(1,pident) %>%
    tidyr::separate(stitle, c("acc", ranks), delim, remove = TRUE)

  if(resolve_ties == "first"){
    top_hit <- top_hit %>%
      dplyr::mutate(row_n = dplyr::row_number()) %>%
      dplyr::top_n(1, row_n) %>% # Break ties by position
      dplyr::select(-row_n) %>%
      dplyr::ungroup()
  } else if(resolve_ties == "all"){
    top_hit <- top_hit %>%
      dplyr::ungroup()
  }
  return(top_hit)
}


# BLAST assign species ----------------------------------------------------


#' Assign species using BLAST
#'
#' @description This is to be used alongside a hierarchial classifier such as IDTAXA or RDP to assign additional species level matches.
#' This is designed to be a more flexible version of dada2's assignSpecies function
#' @param query (Required) Query sequence. Accepts a DNABin object, DNAStringSet object, Character string, or filepath.
#' @param db (Required) Reference sequences to conduct search against. Accepts a DNABin object, DNAStringSet object, Character string, or filepath.
#' If DNAbin, DNAStringSet or character string is provided, a temporary fasta file is used to construct BLAST database
#' @param type (Required) type of search to conduct, default 'blastn'
#' @param identity (Required) Minimum percent identity cutoff. Note that this is calculated using all alignments for each query-subject match.
#' @param coverage (Required) Minimum percent query coverage cutoff. Note that this is calculated using all alignments for each query-subject match.
#' @param evalue (Required) Minimum expect value (E) for saving hits
#' @param max_target_seqs (Required) Number of aligned sequences to keep. Even if you are only looking for 1 top hit keep this higher for calculations to perform properly.
#' @param max_hsp (Required) Maximum number of HSPs (alignments) to keep for any single query-subject pair.
#' @param ranks (Required) The taxonomic ranks contained in the fasta headers
#' @param delim (Required) The delimiter between taxonomic ranks in fasta headers
#' @param args (Optional) Extra arguments passed to BLAST
#' @param quiet (Optional) Whether progress should be printed to console, default is FALSE
#' @param remove_db_gaps Whether gaps should be removed from the fasta file used for the database. Note that makeblastdb can fail if there are too many gaps in the sequence.
#'
#' @return
#' @export
#' @import dplyr
#' @import stringr
#' @importFrom tidyr separate
#'
blast_assign_species <- function(query, db, type="blastn",
                                 identity=97, coverage=95, evalue=1e06, max_target_seqs=5, max_hsp=5,
                                 ranks=c("Kingdom", "Phylum","Class", "Order", "Family", "Genus", "Species"), delim=";",
                                 args=NULL, quiet=FALSE, remove_db_gaps=TRUE){

  #Check input contains species and genus
  if(!any(tolower(ranks) %in% c("species", "genus"))){
    stop("Ranks must include Genus and Species")
  }

  #Conduct BLAST
  top_hit <- blast_top_hit(query = query, db = db, type=type,
                           identity=identity, coverage=coverage, evalue=evalue, max_target_seqs=max_target_seqs, max_hsp=max_hsp,
                           ranks=ranks, delim=delim, resolve_ties="all", args=args, quiet=quiet, remove_db_gaps=remove_db_gaps) %>%
    dplyr::filter(!is.na(Species))

  out <- top_hit %>%
    dplyr::group_by(qseqid) %>%
    dplyr::mutate(spp = Species %>% stringr::str_remove("^.* ")) %>%
    dplyr::reframe(spp = paste(sort(unique(spp)), collapse = "/"), Genus, pident, qcovs, max_score, total_score, evalue) %>%
    dplyr::mutate(binomial = paste(Genus, spp)) %>%
    dplyr::distinct() %>%
    dplyr::add_tally() %>%
    dplyr::mutate(binomial =  dplyr::case_when( #Leave unassigned if conflicted at genus level
      n > 1 ~ as.character(NA),
      n == 1 ~ binomial
    )) %>%
    dplyr::select(OTU = qseqid, Genus, Species = binomial, pident, qcovs, max_score, total_score, evalue)

  return(out)
}
