
# .findExecutable ---------------------------------------------------------


#' Find executable
#'
#' @param exe
#' @param interactive
#'
#' @return
#' @export
#'
#' @examples
.findExecutable <- function(exe, interactive=TRUE) {
  path <- Sys.which(exe)
  if(all(path=="")) {
    if(interactive) stop("Executable for ", paste(exe, collapse=" or "), " not found! Please make sure that the software is correctly installed and, if necessary, path variables are set.", call.=FALSE)
    return(character(0))
  }

  path[which(path!="")[1]]
}


# Install BLAST -----------------------------------------------------------

#' Install BLAST
#'
#' @param url (Optional) Default will search for the latest version
#' URL to retrieve BLAST version from.
#' @param dest.dir (Optional)  Default "bin"
#' Directory to install BLAST within.
#' @force Whether existing installs should be forcefully overwritten
#'
#' @return
#' @export
#'
#' @examples
blast_install <- function(url, dest.dir = "bin", force = FALSE) {

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

  if (!dir.exists(dest.dir)) {
    dir.create(dest.dir) # Create first directory
  }

  blast_version <- basename(url) %>% stringr::str_replace("(-x64)(.*?)(?=$)", "")
  if (dir.exists(paste0(dest.dir, "/", blast_version)) && force == FALSE) {
    message("Skipped as BLAST already exists in directory, to overwrite set force to TRUE")
    return(NULL)
  } else  if (dir.exists(paste0(dest.dir, "/", blast_version)) && force == TRUE) {
    unlink(paste0(dest.dir, "/", blast_version), recursive = TRUE) # Remove old version
  }

  destfile <- file.path(dest.dir, basename(url))
  if (exists(destfile)) {
    file.remove(destfile) # Remove old zip file
  }

  #Download and unzip
  httr::GET(url, httr::write_disk(destfile, overwrite=TRUE))
  utils::untar(destfile, exdir = dest.dir)
  file.remove(destfile)

  #Set new $Paths variable for mac & linux
  if(localos == "Darwin" | localos == "unix"){
    old_path <- Sys.getenv("PATH")
    install_path <- list.dirs(dest.dir, full.names = TRUE)[str_detect(list.dirs(dest.dir, full.names = TRUE),"/bin$")]
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
#'
#' @return
#' @export
#'
#' @examples
makeblastdb <- function (file, dbtype = "nucl", args = NULL, quiet = FALSE) {
  time <- Sys.time() # get time
  .findExecutable("makeblastdb") # Check blast is installed
  if (is.null(args)){args <- ""}
  if (stringr::str_detect(file, ".gz")) {
    message("Unzipping file")
    compressed <- TRUE
    R.utils::gunzip(file, remove=FALSE)
    file <- stringr::str_replace(file, ".gz", "")
  }else (compressed <- FALSE)
  results <- system2(command = .findExecutable("makeblastdb"),
                     args = c("-in", file, "-dbtype", dbtype, args),
                     wait = TRUE,
                     stdout = TRUE)
  time <- Sys.time() - time
  if (compressed) {file.remove(file)}
  if (!quiet) (message(paste0("made BLAST DB in ", format(time, digits = 2))))

}


#' Show BLAST parameters
#'
#' @param type (Required) Which BLAST function to display help page for
#'
#' @return
#' @export
#'
#' @examples
blast_params <- function(type = "blastn") {
  system(paste(.findExecutable(c(type)), "-help"))
}


# BLAST -------------------------------------------------------------------

#' Run BLAST search
#'
#' @param query (Required) Query sequence. Accepts a DNABin object, DNAStringSet object, Character string, or filepath.
#' @param db (Required) Reference sequences to conduct search against. Accepts a DNABin object, DNAStringSet object, Character string, or filepath.
#' If DNAbin, DNAStringSet or character string is provided, a temporary fasta file is used to construct BLAST database
#' @param type (Required) type of search to conduct, default 'blastn'
#' @param evalue (Required) Minimum evalue from search
#' @param args (Optional) Extra arguments passed to BLAST
#' @param quiet (Optional) Whether progress should be printed to console, default is FALSE
#'
#' @return
#' @export
#'
#' @examples
blast <- function (query, db, type="blastn", evalue = 1e-6, output_format = "tabular", args=NULL, quiet=FALSE){
  time <- Sys.time() # get time
  # Create temp files
  tmp <- tempdir()
  tmpquery <- paste0(tmp, "/tmpquery.fa")
  tmpdb <- paste0(tmp, "/tmpquery.fa")

  .findExecutable(type) # check blast is installed

  # Check outfmt
  if(output_format=="tabular"){
    outfmt <- 6
  } else if (output_format=="lulu"){
    outfmt <- "'6 qseqid sseqid pident'"
  } else if(!output_format==6 && is.numeric(output_format)){
    outfmt <- output_format
  } else if(!output_format %in% c("tabular", "lulu") && is.character(output_format)){
    outfmt <- paste0("'",output_format, "'")
  }

  # Database
  if(inherits(db, "DNAbin")){
    if (!quiet) { message("Database input is DNAbin: Creating temporary blast database") }
    insect::writeFASTA(db, tmpdb)
    makeblastdb(tmpdb)
    db <- tmpdb
  } else if (inherits(db, "DNAString") | inherits(db, "DNAStringSet")){
    if (!quiet) { message("Database input is DNAStringSet: Creating temporary blast database") }
    Biostrings::writeXStringSet(db, tmpdb)
    makeblastdb(tmpdb)
    db <- tmpdb
  } else if (inherits(db, "character") &&  all(stringr::str_to_upper(stringr::str_split(db,"")[[1]]) %in% Biostrings::DNA_ALPHABET)) { # Handle text input
    if (!quiet) { message("Database input is character string: Creating temporary blast database") }
    if (nchar(db[1]) == 1) {db <- paste0(db, collapse = "")}
    db <- insect::char2dna(db)
    insect::writeFASTA(db, tmpdb)
    makeblastdb(tmpdb)
    db <- tmpdb
  } else if (inherits(db, "character") &&  file.exists(file.path(db))){ # Handle filename
    makeblastdb(db)
    db <- stringr::str_replace(db, ".gz", "")
    db <- db
  } else {
    stop("Invalid BLAST database")
  }

  # Query
  if(inherits(query, "DNAbin")){
    if (!quiet) { message("Query is DNAbin: Creating temporary fasta file") }
    insect::writeFASTA(query, tmpquery)
    input <- tmpquery
  } else if (inherits(query, "DNAString") | inherits(query, "DNAStringSet")){
    if (!quiet) { message("Query is DNAString: Creating temporary fasta file") }
    Biostrings::writeXStringSet(query, tmpquery)
    input <- tmpquery
  }else if (inherits(query, "character") &&  all(str_to_upper(str_split(query,"")[[1]]) %in% Biostrings::DNA_ALPHABET)) { # Handle text query
    if (!quiet) { message("Query is character string: Creating temporary fasta file") }
    if (nchar(query[1]) == 1) {query <- paste0(query, collapse = "")}
    query <- insect::char2dna(query)
    insect::writeFASTA(query, tmpquery)
    input <- tmpquery
  } else if (inherits(query, "character") &&  file.exists(file.path(query))){ # Handle filenames
    input <- query
  }else {
    stop("Invalid BLAST query")
  }

  #  Conduct BLAST search
  if (!quiet) { message("Running BLAST") }
  results <- system2(command = .findExecutable(type),
                     args = c("-db", db,
                              "-query", input,
                              "-outfmt ", outfmt,
                              "-evalue", evalue,
                              "-ungapped", args),
                     wait = TRUE,
                     stdout = TRUE)

  # Parse BLAST results
  if(output_format=="tabular"){
    out <- results %>%
      tibble::enframe() %>%
      tidyr::separate(col = value,
                      into = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"),
                      sep = "\t",
                      convert = TRUE)
  } else if(output_format =="lulu"){
    out <- results %>%
      tibble::enframe() %>%
      tidyr::separate(col = value,
                      into = c("qseqid", "sseqid", "pident"),
                      sep = "\t",
                      convert = TRUE)
  } else{
    message("Warning, result parsing is currently only supported for output_format = 'tabular'
            and output_format = 'lulu', returning raw results")
    out <- results %>%
      tibble::enframe()
  }
  time <- Sys.time() - time
  if (!quiet) (message(paste0("finished BLAST in ", format(time, digits = 2))))

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
#' @param threshold (Required) Minimum identity threshold to accept
#' @param taxranks (Required) The taxonomic ranks contained in the fasta headers
#' @param delim (Required) The delimiter between taxonomic ranks in fasta headers
#' @param args (Optional) Extra arguments passed to BLAST
#' @param quiet (Optional) Whether progress should be printed to console, default is FALSE
#'
#' @return
#' @export
#'
#' @examples
blast_top_hit <- function(query, db, type="blastn", threshold=90, taxranks=c("Kingdom", "Phylum","Class", "Order", "Family", "Genus", "Species"), delim=";", args="-max_target_seqs 5", quiet=FALSE ){
  #Conduct BLAST
  result <- blast(query=query, type=type, db=db, args=args)
  #Subset to top hit
  top_hit <- result %>%
    filter(pident > threshold) %>%
    group_by(qseqid) %>%
    top_n(1, bitscore) %>%
    top_n(1, row_number(name)) %>% # Break ties by position
    separate(sseqid, c("acc",taxranks), delim)

}
