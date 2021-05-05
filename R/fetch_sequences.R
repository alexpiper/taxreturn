# Genbank functions -----------------------------------------------
# Some useful entrez queries
#' all `[filter]` 	Retrieves everthing
#' Specified `[property]` 	Formal binomial and trinomial
#' at or below species level `[property]`
#' family `[rank]` 	Rank-based query
#' taxonomy genome `[filter]` 	Taxa with a direct link to a genome sequence
#' 2009/10/21:2020 `[date]` 	Date-bounded query
#' mammalia `[subtree]` 	All taxa within the Mammalia
#' extinct `[property]` 	Extinct organisms
#' Terminal `[property]` 	Terminal nodes in the tree
#' loprovencyclife `[filter]` 	Entries with LinkOut links to the Encyclopedia of Life
#'

#' Genbank search function
#'
#' @param x A taxon name or vector of taxa to download sequences for
#' @param database The database to download from. For NCBI GenBank this currently onlt accepts the arguments 'nuccore' or 'genbank' which is an alias for nuccore.
#' @param marker The barcode marker used as a search term for the database.
#' If this is set to "mitochondria" or "mitochondrion" it will download full mitochondrial genomes. If set to "genome" it will download entire genomes only.
#' @param quiet Whether progress should be printed to the console.
#' @param output The output format for the taxonomy in fasta headers.
#' Options include "h" for full heirarchial taxonomy (SeqID;Domain;Phylum;Class;Order;Family;Genus;Species),
#' "binom" for just genus species binomials (SeqID;Genus Species),
#' "bold" for BOLD taxonomic ID only (SeqID;BoldTaxID),
#' "gb" for genbank taxonomic ID (SeqID;GBTaxID),
#' "gb-binom" which outputs Genus species binomials, as well as genbank taxonomic ID's, and translates all BOLD taxonomic ID's to genbank taxonomic ID's in the process
#' or "standard" which outputs the default format for each database. For genbank this is `Accession Sequence definition`
#' @param min_length The minimum length of sequences to download
#' @param max_length The maximum length of sequences to download
#' @param chunk_size Split up the query into chunks of this size to avoid overloading API servers.
#' if left NULL, the default will be 10,000 for regular queries, 1,000 if marker is "mitochondria", and 1 if marker is "genome"
#' @param compress  Option to compress output fasta files using gzip
#' @param force Option ot overwright files if they already exist
#' @param out_dir Output directory to write fasta files to
#'
#'
#' @import dplyr
#' @import stringr
#' @import purrr
#' @importFrom biofiles gbRecord
#' @importFrom biofiles getAccession
#' @importFrom biofiles getDefinition
#' @importFrom biofiles getTaxonomy
#' @importFrom biofiles getOrganism
#' @importFrom biofiles getSequence
#' @importFrom biofiles dbxref
#' @importFrom rentrez entrez_search
#' @importFrom rentrez entrez_fetch
#' @importFrom tidyr unite
#' @importFrom Biostrings fasta.index
#' @importFrom Biostrings writeXStringSet
#'
#' @return
#' @details
#'
search_genbank <- function(x, database = "nuccore", marker = c("COI[GENE]", "CO1[GENE]", "COX1[GENE]"), quiet = FALSE, output = "h",
                     min_length = 1, max_length = 2000, chunk_size=NULL, out_dir = NULL,
                     compress = FALSE, force=FALSE) {

  # function setup
  time <- Sys.time() # get time

  if(!database %in% c("nuccore", "genbank")){
    stop("database is invalid: only nuccore and genbank is currently supported")
  }

  if (!output %in% c("standard", "h", "binom", "gb", "bold", "gb-binom")) {
    stop(paste0(output, " has to be one of: 'standard', 'h','binom','bold', 'gb' or 'gb-binom', see help page for more details"))
  }
  if (!any(tolower(marker) %in% c("mitochondria", "mitochondrion", "genome"))) {
    marker <- paste(marker, collapse=" OR ")
  }
  #Define directories
  if (is.null(out_dir)) {
    out_dir <- database
    if (!quiet) (message(paste0("No input out_dir given, saving output file to: ", out_dir)))
  }
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }

  # Create output file
  name <- marker %>%
     stringr::str_remove_all(pattern="\\[GENE]") %>%
     stringr::str_remove_all(pattern="OR ") %>%
     stringr::str_replace_all(pattern=" ", replacement ="_")

  if (compress) {
    out_file <- paste0(normalizePath(out_dir), "/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", name, ".fa.gz")
  } else if (!compress) {
    out_file <- paste0(normalizePath(out_dir), "/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", name, ".fa")
  }

  #Check if output exists
  if (file.exists(out_file) && force==TRUE) {
    file.remove(out_file)
    cat("", file=out_file)
  } else if (file.exists(out_file) && force==FALSE){
    stop(paste0(out_file, " exists, set force = TRUE to overwrite"))
  }

  if(database=="genbank"){
    database="nuccore"
  }

  # Main function
  tryCatch(
    {
      if (!tolower(marker) %in% c("mitochondria", "mitochondrion", "genome")) {
        searchQ <- paste("(", x, "[ORGN])", " AND (", paste(c(marker), collapse = " OR "), ") AND ", min_length, ":", max_length, "[Sequence Length]", sep = "")
        if(is.null(chunk_size)) {chunk_size=10000}
      } else if (tolower(marker) %in% c("mitochondria", "mitochondrion")) {
        message(paste0("Input marker is ", marker, ", Downloading full mitochondrial genomes"))
        searchQ <- paste("(", x, "[ORGN])", " AND mitochondrion[filter] AND genome", sep = "")
        if(is.null(chunk_size)) {chunk_size=1000}
      } else if (tolower(marker) %in% c("genome")){
        message(paste0("Input marker is ", marker, ", Downloading full genomes"))
        searchQ <- paste("(", x, "[ORGN])", " AND complete genome[title]", sep = "")
        if(is.null(chunk_size)) {chunk_size=1}
      }
      if(!quiet){message("Searching genbank with query:", searchQ)}

      search_results <- rentrez::entrez_search(db = database, term = searchQ, retmax = 9999999, use_history = TRUE)

      if (search_results$count > 0 & !is.na(search_results$count)) {
        if (!quiet) {message(paste0(search_results$count, " Sequences to be downloaded for: ", searchQ))}

        l <- 1
        start <- 0

        chunks <- length(search_results$ids) / chunk_size
        if (!is.integer(chunks)) {
          chunks <- as.integer(length(search_results$ids) / chunk_size) + 1
        }

        for (l in 1:chunks) {

          # Fetch gb flat files
          dl <- rentrez::entrez_fetch(db = database, web_history = search_results$web_history, rettype = "gb", retmax = chunk_size, retstart = start)
          gb <- biofiles::gbRecord(rcd = textConnection(dl))

          # Hierarchial output
          if (output == "standard") {
            names <- paste0(biofiles::getAccession(gb), " ", biofiles::getDefinition(gb))
          } else if (output == "h") {
            lineage <- biofiles::getTaxonomy(gb) %>%
              stringr::str_remove(".$") %>%
              stringr::str_split_fixed(pattern = ";", n = Inf) %>%
              trimws(which = "both") %>%
              as.data.frame() %>%
              dplyr::mutate(Species = biofiles::getOrganism(gb)) %>%
              tidyr::unite("names", everything(), sep = ";") %>%
              dplyr::mutate(names = names %>%
                              stringr::str_replace_all(pattern = " ", replacement = "_")%>%
                              stringr::str_replace_all(pattern = ";;", replacement = ";"))
            names <- paste0(biofiles::getAccession(gb), ";", lineage$names)

            # Genbank taxID output
          } else if (output == "gb") {
            names <- cat_acctax(gb)
          } else if (output == "binom") {
            names <- paste0(biofiles::getAccession(gb), ";", biofiles::getOrganism(gb))%>%
              stringr::str_replace_all(pattern = " ", replacement = "_")
          } else if (output == "gb-binom") {
            names <- paste0(cat_acctax(gb), ";", biofiles::getOrganism(gb))
          }
          # Output FASTA
          seqs <- biofiles::getSequence(gb) # Something failing here - getting much less than the proper amount
          if(all(is.na(names(seqs)))) {
            names(seqs) <- biofiles::getAccession(gb)
          }
          #Check if names match
          names(seqs) <- names[names %>%
                                 stringr::str_remove(pattern=";.*$") %>%
                                 stringr::str_remove(pattern="\\|.*$") %>%
                                 stringr::str_remove(" .*$")
                               %in% names(seqs)]
          if (compress == TRUE) {
            Biostrings::writeXStringSet(seqs, out_file, format = "fasta", compress = "gzip", width = 20000, append = TRUE)
          } else if (compress == FALSE) {
            Biostrings::writeXStringSet(seqs, out_file, format = "fasta", width = 20000, append = TRUE)
          }

          if (!quiet) (message("Chunk", l, " of ", chunks, " downloaded\r"))
          start <- start + chunk_size
          Sys.sleep(2.5)
          if (l >= chunks) {
            time <- Sys.time() - time
            if (!quiet) (message("All chunks downloaded \r"))
          }
        }
      }
    },
    error = function(e) NULL
  )

  #Count number of downloaded sequences
  if(file.exists(out_file)){
    counter <- nrow(Biostrings::fasta.index(out_file))
  } else {
    counter <- 0
  }

  # Count total sequences that should have been downloaded
  if(exists("search_results") && length(search_results$ids) > 1 ){
    total_counter <- length(search_results$ids)
  } else {
    total_counter <- 0
  }

  res <- data.frame(
    taxon = x,
    seqs_total = total_counter,
    seqs_downloaded = counter,
    marker = marker,
    database = database,
    time = format(time, digits = 2)
  )
  return(res)
}

#' Genbank subsampling
#'
#' @param x A taxon name or vector of taxa to download sequences for
#' @param database The database to download from. For NCBI GenBank this currently onlt accepts the arguments 'nuccore' or 'genbank' which is an alias for nuccore.
#' @param marker The barcode marker used as a search term for the database
#' @param quiet Whether progress should be printed to the console.
#' @param output The output format for the taxonomy in fasta headers.
#' Options include "h" for full heirarchial taxonomy (SeqID;Domain;Phylum;Class;Order;Family;Genus;Species),
#' "binom" for just genus species binomials (SeqID;Genus Species),
#' "bold" for BOLD taxonomic ID only (SeqID;BoldTaxID),
#' "gb" for genbank taxonomic ID (SeqID;GBTaxID),
#' "gb-binom" which outputs Genus species binomials, as well as genbank taxonomic ID's, and translates all BOLD taxonomic ID's to genbank taxonomic ID's in the process
#' or "standard" which outputs the default format for each database. For genbank this is `Accession Sequence definition`
#' @param min_length The minimum length of sequences to download
#' @param max_length The maximum length of sequences to download
#' @param subsample (Numeric) return a random subsample of sequences from the search.
#' @param chunk_size Split up the query into chunks of this size to avoid overloading API servers.
#' if left NULL, the default will be 10,000 for regular queries, 1,000 if marker is "mitochondria", and 1 if marker is "genome"
#' @param compress  Option to compress output fasta files using gzip
#' @param force Option ot overwright files if they already exist
#' @param out_dir Output directory to write fasta files to
#'
#'
#' @import dplyr
#' @import stringr
#' @import purrr
#' @importFrom biofiles gbRecord
#' @importFrom biofiles getAccession
#' @importFrom biofiles getDefinition
#' @importFrom biofiles getTaxonomy
#' @importFrom biofiles getOrganism
#' @importFrom biofiles getSequence
#' @importFrom biofiles dbxref
#' @importFrom rentrez entrez_search
#' @importFrom rentrez entrez_post
#' @importFrom rentrez entrez_fetch
#' @importFrom tidyr unite
#' @importFrom Biostrings fasta.index
#' @importFrom Biostrings writeXStringSet
#'
#' @return
#'
#'
search_genbank_subsample <- function(x, database = "nuccore", marker = c("COI[GENE]", "CO1[GENE]", "COX1[GENE]"),
                               quiet = FALSE, output = "h", min_length = 1, max_length = 2000,
                               subsample=1000, chunk_size=300, compress = FALSE, force=FALSE, out_dir = NULL) {
  # function setup
  time <- Sys.time() # get time

  if(!database %in% c("nuccore", "genbank")){
    stop("database is invalid: only nuccore and genbank is currently supported")
  }

  if (!output %in% c("standard", "h", "binom", "gb", "bold", "gb-binom")) {
    stop(paste0(output, " has to be one of: 'standard', 'h','binom','bold', 'gb' or 'gb-binom', see help page for more details"))
  }
  if (!any(tolower(marker) %in% c("mitochondria", "mitochondrion", "genome"))) {
    marker <- paste(marker, collapse=" OR ")
  }
  #Define directories
  if (is.null(out_dir)) {
    out_dir <- database
    if (!quiet) (message(paste0("No input out_dir given, saving output file to: ", out_dir)))
  }
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }

  # Create output file
  name <- marker %>%
    stringr::str_remove_all(pattern="\\[GENE]") %>%
    stringr::str_remove_all(pattern="OR ") %>%
    stringr::str_replace_all(pattern=" ", replacement ="_")

  if (compress) {
    out_file <- paste0(normalizePath(out_dir), "/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", name, ".fa.gz")
  } else if (!compress) {
    out_file <- paste0(normalizePath(out_dir), "/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", name, ".fa")
  }

  #Check if output exists
  if (file.exists(out_file) && force==TRUE) {
    file.remove(out_file)
    cat("", file=out_file)
  } else if (file.exists(out_file) && force==FALSE){
    stop(paste0(out_file, " exists, set force = TRUE to overwrite"))
  }

  if(database=="genbank"){
    database="nuccore"
  }

  # Main function
  tryCatch(
    {
      # Genbank Search
      if (!tolower(marker) %in%c("mitochondria", "mitochondrion", "genome")) {
        searchQ <- paste("(", x, "[ORGN])", " AND (", paste(c(marker), collapse = " OR "), ") AND ", min_length, ":", max_length, "[Sequence Length]", sep = "")
        chunk_size=5000
      } else if (tolower(marker) %in% c("mitochondria", "mitochondrion")) {
        message(paste0("Input marker is ", marker, ", Downloading full mitochondrial genomes"))
        searchQ <- paste("(", x, "[ORGN])", " AND mitochondrion[filter] AND genome", sep = "")
        chunk_size=500
      } else if (tolower(marker) %in% c("genome", "mitochondrion")){
        message(paste0("Input marker is ", marker, ", Downloading full genomes"))
        searchQ <- paste("(", x, "[ORGN])", " AND complete genome[title]", sep = "")
        chunk_size=1
      }
      if(!quiet){message("Searching genbank with query:", searchQ)}

      search_results <- rentrez::entrez_search(db = database, term = searchQ, retmax = 9999999, use_history = TRUE)

      if (search_results$count > 0 & !is.na(search_results$count)) {

        if (!quiet) {message(paste0(subsample, " of: ", search_results$count, " sequences to be downloaded for: ", searchQ))}

        ids <- sample(search_results$ids, subsample )
        chunks <- split(ids, ceiling(seq_along(ids)/chunk_size))

        l <- 1

        for (l in 1:length(chunks)) {

          #Upload ids
          upload <- rentrez::entrez_post(db=database, id=chunks[[l]])

          # Fetch gb flat files
          dl <- rentrez::entrez_fetch(db = database, web_history = upload, rettype = "gb", retmax = 10000)
          gb <- biofiles::gbRecord(rcd = textConnection(dl))

          # Hierarchial output
          if (output == "standard") {
            names <- paste0(biofiles::getAccession(gb), " ", biofiles::getDefinition(gb))
          } else if (output == "h") {
            lineage <- biofiles::getTaxonomy(gb) %>%
              str_split_fixed(pattern = ";", n = Inf) %>%
              trimws(which = "both") %>%
              as_tibble() %>%
              dplyr::mutate(Species = biofiles::getOrganism(gb)) %>%
              dplyr::mutate(Genus = str_replace(V15, pattern = "[.]", replacement = "")) %>%
              tidyr::unite("names", c("V1", "V2", "V4", "V6", "V10", "V14", "Genus", "Species"), sep = ";") %>%
              dplyr::mutate(names = str_replace(names, pattern = " ", replacement = "_"))
            names <- paste0(names(biofiles::getSequence(gb)), ";", lineage$names)

            # Genbank taxID output
          } else if (output == "gb") {
            names <- cat_acctax(gb)
          } else if (output == "binom") {
            names <- paste0(biofiles::getAccession(gb), ";", biofiles::getOrganism(gb))
          } else if (output == "gb-binom") {
            names <- paste0(cat_acctax(gb), ";", biofiles::getOrganism(gb))
          }
          #output seqs
          seqs <- biofiles::getSequence(gb) # Something failing here - getting much less than the proper amount
          if(all(is.na(names(seqs)))) {
            names(seqs) <- biofiles::getAccession(gb)
          }
          #Check if names match
          names(seqs) <- names[names %>%
                                 stringr::str_remove(pattern=";.*$") %>%
                                 stringr::str_remove(pattern="\\|.*$") %>%
                                 stringr::str_remove(" .*$")
                               %in% names(seqs)]
          if (compress == TRUE) {
            Biostrings::writeXStringSet(seqs, out_file, format = "fasta", compress = "gzip", width = 20000, append = TRUE)
          } else if (compress == FALSE) {
            Biostrings::writeXStringSet(seqs, out_file, format = "fasta", width = 20000, append = TRUE)
          }

          if (!quiet) (message("Chunk", l, " of ", length(chunks), " downloaded\r"))
          Sys.sleep(2.5)
          if (l >= length(chunks)) {
            time <- Sys.time() - time
            if (!quiet) (message("All chunks downloaded \r"))
          }
        }
      }
    },
    error = function(e) NULL
  )

  #Count number of downloaded sequences
  if(file.exists(out_file)){
    counter <- nrow(Biostrings::fasta.index(out_file))
  } else {
    counter <- 0
  }

  # Count total sequences that should have been downloaded
  if(exists("search_results") && length(search_results$ids) > 1 ){
    total_counter <- length(search_results$ids)
  } else {
    total_counter <- 0
  }

  res <- data.frame(
    taxon = x,
    seqs_total = total_counter,
    seqs_downloaded = counter,
    marker = marker,
    database = database,
    time = format(time, digits = 2)
  )
  return(res)
}


#' GbUpdate
#'
#' @param x A taxon name or vector of taxa to download sequences for
#' @param fasta A fasta file or list of fasta files to check for existing sequence accessions
#' @param database The database to download from. For NCBI GenBank this currently onlt accepts the arguments 'nuccore' or 'genbank' which is an alias for nuccore.
#' @param marker The barcode marker used as a search term for the database. If this is set to "mitochondria" it will download full mitochondrial genomes.
#' @param quiet Whether progress should be printed to the console.
#' @param output The output format for the taxonomy in fasta headers.
#' Options include "h" for full heirarchial taxonomy (SeqID;Domain;Phylum;Class;Order;Family;Genus;Species),
#' "binom" for just genus species binomials (SeqID;Genus Species),
#' "bold" for BOLD taxonomic ID only (SeqID;BoldTaxID),
#' "gb" for genbank taxonomic ID (SeqID;GBTaxID),
#' or "gb-binom" which outputs Genus species binomials, as well as genbank taxonomic ID's, and translates all BOLD taxonomic ID's to genbank taxonomic ID's in the process
#' @param suffix The suffix to add to newly downloaded files. Defaults to 'updates'
#' @param min_length The minimum length of sequences to download
#' @param max_length The maximum length of sequences to download
#' @param chunk_size The size of the chunked searches to conduct.
#' Warning, chunk sizes over 300 can be too big for the NCBI servers.
#' @param compress  Option to compress output fasta files using gzip
#' @param force Option ot overwright files if they already exist
#' @param out_dir Output directory to write fasta files to
#' @param multithread Whether multithreading should be used, if TRUE the number of cores will be automatically detected, or provided a numeric vector to manually set the number of cores to use
#' @param progress Whether a progress bar should be shown. This reduces speed and is false by default.
#'
#' @return
#' @import stringr
#' @import future
#' @import furrr
#' @importFrom biofiles gbRecord
#' @importFrom biofiles getAccession
#' @importFrom biofiles getDefinition
#' @importFrom biofiles getTaxonomy
#' @importFrom biofiles getOrganism
#' @importFrom biofiles getSequence
#' @importFrom biofiles dbxref
#' @importFrom rentrez entrez_search
#' @importFrom rentrez entrez_post
#' @importFrom rentrez entrez_fetch
#' @importFrom Biostrings fasta.index
#' @importFrom Biostrings writeXStringSet
#'
update_genbank <- function(x, fasta, database = "nuccore", marker = c("COI[GENE]", "CO1[GENE]", "COX1[GENE]"), quiet = FALSE, output = "h", suffix="updates",
                     min_length = 1, max_length = 2000, chunk_size=300, out_dir = NULL,
                     compress = FALSE, force=FALSE, multithread = FALSE, progress=FALSE){

  # function setup
  time <- Sys.time() # get time

  if(!database %in% c("nuccore", "genbank")){
    stop("database is invalid: only nuccore and genbank is currently supported")
  }

  #Check if all inputs exists
  if (!any(file.exists(fasta))) {
    stop("Not all fasta files provided can be found, check the diferectory you have provided")
  }

  if (!output %in% c("standard", "h", "binom", "gb", "bold", "gb-binom")) {
    stop(paste0(output, " has to be one of: 'standard', 'h','binom','bold', 'gb' or 'gb-binom', see help page for more details"))
  }
  if (!any(tolower(marker) %in% c("mitochondria", "mitochondrion", "genome"))) {
    marker <- paste(marker, collapse=" OR ")
  }
  #Define directories
  if (is.null(out_dir)) {
    out_dir <- database
    if (!quiet) (message(paste0("No input out_dir given, saving output file to: ", out_dir)))
  }
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }

  # Create output file
  name <- marker %>%
    stringr::str_remove_all(pattern="\\[GENE]") %>%
    stringr::str_remove_all(pattern="OR ") %>%
    stringr::str_replace_all(pattern=" ", replacement ="_")

  if (compress) {
    out_file <- paste0(normalizePath(out_dir), "/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", name, "_", suffix, ".fa.gz")
  } else if (!compress) {
    out_file <- paste0(normalizePath(out_dir), "/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", name, "_", suffix, ".fa")
  }

  #Check if output exists
  if (file.exists(out_file) && force==TRUE) {
    file.remove(out_file)
    cat("", file=out_file)
  } else if (file.exists(out_file) && force==FALSE){
    stop(paste0(out_file, " exists, set force = TRUE to overwrite"))
  }

  if(database=="genbank"){
    database="nuccore"
  }

  # Get accessions from existing fastas
  current <- acc_from_fasta(fasta)
  if(!quiet){message(length(current), " unique accessions in fastas")}

  # Genbank Search
  if (!tolower(marker) %in%c("mitochondria", "mitochondrion", "genome")) {
    searchQ <- paste("(", x, "[ORGN])", " AND (", paste(c(marker), collapse = " OR "), ") AND ", min_length, ":", max_length, "[Sequence Length]", sep = "")
    if(is.null(chunk_size)) {chunk_size=10000}
  } else if (tolower(marker) %in% c("mitochondria", "mitochondrion")) {
    message(paste0("Input marker is ", marker, ", Downloading full mitochondrial genomes"))
    searchQ <- paste("(", x, "[ORGN])", " AND mitochondrion[filter] AND genome", sep = "")
    if(is.null(chunk_size)) {chunk_size=1000}
  } else if (tolower(marker) %in% c("genome", "mitochondrion")){
    message(paste0("Input marker is ", marker, ", Downloading full genomes"))
    searchQ <- paste("(", x, "[ORGN])", " AND complete genome[title]", sep = "")
    if(is.null(chunk_size)) {chunk_size=1}
  }
  if(!quiet){message("Searching genbank with query:", searchQ)}

  search_results <- rentrez::entrez_search(db = database, term = searchQ, retmax = 9999999, use_history = TRUE)
  ids <- search_results$ids

  #Convert search result GIDs to accession
  if(length(ids) > 0 ){
    accs <- gid_to_acc(ids, database = database, chunk_size = chunk_size,  multithread=multithread, progress = progress) %>%
    stringr::str_remove(".[0-9]$")
  } else{
    stop("search returned no hits")
  }
  # Find any missing accessions
  newsearch <- setdiff(accs, current)

  # Download missing accessions

  if (length(newsearch) > 0) {

    if (!quiet) {message(paste0(length(newsearch), " sequences to be downloaded for: ", searchQ))}

    chunks <- split(newsearch, ceiling(seq_along(newsearch)/chunk_size))

    # setup multithreading
    setup_multithread(multithread = multithread, quiet=quiet)

    #Main function
    furrr::future_map(chunks, function(x){
      upload <- rentrez::entrez_post(db=database, id=x)
      dl <- rentrez::entrez_fetch(db = database, web_history = upload, rettype = "gb", retmax = chunk_size)
      gb <- biofiles::gbRecord(rcd = textConnection(dl))

      # Hierarchial output
      if (output == "standard") {
        names <- paste0(biofiles::getAccession(gb), " ", biofiles::getDefinition(gb))
      } else if (output == "h") {
        lineage <- biofiles::getTaxonomy(gb) %>%
          str_split_fixed(pattern = ";", n = Inf) %>%
          trimws(which = "both") %>%
          as_tibble() %>%
          dplyr::mutate(Species = biofiles::getOrganism(gb)) %>%
          dplyr::mutate(Genus = str_replace(V15, pattern = "[.]", replacement = "")) %>%
          tidyr::unite("names", c("V1", "V2", "V4", "V6", "V10", "V14", "Genus", "Species"), sep = ";") %>%
          dplyr::mutate(names = str_replace(names, pattern = " ", replacement = "_"))
        names <- paste0(names(biofiles::getSequence(gb)), ";", lineage$names)

        # Genbank taxID output
      } else if (output == "gb") {
        names <- cat_acctax(gb)
      } else if (output == "binom") {
        names <- paste0(biofiles::getAccession(gb), ";", biofiles::getOrganism(gb))
      } else if (output == "gb-binom") {
        names <- paste0(cat_acctax(gb), ";", biofiles::getOrganism(gb))
      }
      seqs <- biofiles::getSequence(gb) # Something failing here - getting much less than the proper amount
      if(all(is.na(names(seqs)))) {
        names(seqs) <- biofiles::getAccession(gb)
      }
      #Check if names match
      names(seqs) <- names[names %>%
                             stringr::str_remove(pattern=";.*$") %>%
                             stringr::str_remove(pattern="\\|.*$") %>%
                             stringr::str_remove(" .*$")
                           %in% names(seqs)]
      if (compress == TRUE) {
        Biostrings::writeXStringSet(seqs, out_file, format = "fasta", compress = "gzip", width = 20000, append = TRUE)
      } else if (compress == FALSE) {
        Biostrings::writeXStringSet(seqs, out_file, format = "fasta", width = 20000, append = TRUE)
      }
      invisible(NULL)
    }, .progress = progress)

  }
  #Close all workers
  future::plan(future::sequential)

  #Count number of downloaded sequences
  if(file.exists(out_file)){
    counter <- nrow(Biostrings::fasta.index(out_file))
  } else {
    counter <- 0
  }

  # Count total sequences that should have been downloaded
  if(exists("search_results") && length(search_results$ids) > 1 ){
    total_counter <- length(search_results$ids)
  } else {
    total_counter <- 0
  }

  res <- data.frame(
    taxon = x,
    seqs_total = total_counter,
    seqs_downloaded = counter,
    marker = marker,
    database = database,
    time = format(time, digits = 2)
  )
  return(res)
}

#' Helper function for processing taxids
#'
#' @param x a gb object
#' @importFrom biofiles getAccession
#' @importFrom biofiles dbxref
cat_acctax <- function(x) {
  if(length(biofiles::getAccession(x)) > 1){
    taxid <- purrr::map(x, biofiles::dbxref, "taxon")
    tax_chr <- purrr::map_chr(taxid, ~{
      as.character(.x[1, 1])})
    attributes(tax_chr) <- NULL
    taxout <- paste0(attributes(taxid)$names, "|", tax_chr)
  } else if(length(biofiles::getAccession(x)) == 1){
    taxout <- paste0(biofiles::getAccession(x), "|", as.character(biofiles::dbxref(x[1], "taxon")))
  }
  return(taxout)
}


# BOLD functions --------------------------------------------------------

#' seach_bold
#' @param x A taxon name or vector of taxa to download sequences for
#' @param marker the barcode marker used as a search term for the database
#' @param quiet Whether progress should be printed to the console.
#' @param output the output format for the taxonomy in fasta headers.
#' Options include "h" for full heirarchial taxonomy (SeqID;Domain;Phylum;Class;Order;Family;Genus;Species),
#' "binom" for just genus species binomials (SeqID;Genus Species),
#' "bold" for BOLD taxonomic ID only (SeqID;BoldTaxID),
#' "gb" for genbank taxonomic ID (SeqID;GBTaxID),
#' "gb-binom" which outputs Genus species binomials, as well as genbank taxonomic ID's, and translates all BOLD taxonomic ID's to genbank taxonomic ID's in the process
#' or "standard" which outputs the default format for each database. For bold this is `sampleid|species name|markercode|genbankid`
#' @param out_file The file to write to, if empty it defaults to the search term
#' @param compress Option to compress output fasta files using gzip
#' @param force Option ot overwright files if they already exist
#' @param out_dir Output directory to write fasta files to
#' @param db (Optional) a database file generated using `taxreturn::get_ncbi_taxonomy()` or `taxreturn::get_ott_taxonomy()`
#'
#' @import dplyr
#' @import stringr
#' @importFrom bold bold_seqspec
#' @importFrom tidyr unite
#' @importFrom tidyr drop_na
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings fasta.index
#' @importFrom Biostrings writeXStringSet
#' @importFrom methods is
#'
#' @return
#'
#'
seach_bold <- function(x, marker = "COI-5P", quiet = FALSE, output = "h",
                       out_file = NULL, compress = FALSE, force=FALSE,
                       out_dir = NULL, db=NULL) {

  # function setup
  time <- Sys.time() # get time

  if (!output %in% c("standard", "h", "binom", "gb", "bold", "gb-binom")) {
    stop(paste0(output, " has to be one of: 'standard', 'h','binom','bold', 'gb' or 'gb-binom', see help page for more details"))
  }

  # Get NCBI taxonomy database if NCBI format outputs are desired
  if (is.null(db) && output %in% c("gb", "gb-binom")) {
    db <- get_ncbi_taxonomy(include_synonyms = TRUE, force=FALSE)
  }

  #Define directories
  if (is.null(out_dir)) {
    out_dir <- "bold"
    if (!quiet) (message(paste0("No input out_dir given, saving output file to: ", out_dir)))
  }
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }

  # Create output files
  if (compress) {
    out_file <- paste0(normalizePath(out_dir), "/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", marker, ".fa.gz")
  } else if (!compress) {
    out_file <- paste0(normalizePath(out_dir), "/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", marker, ".fa")
  }

  #Check if output exists
  if (file.exists(out_file) && force==TRUE) {
    file.remove(out_file)
    cat("", file=out_file)
  } else if (file.exists(out_file) && force==FALSE){
    stop(paste0(out_file, " exists, set force = TRUE to overwrite"))
  }

  # Bold search
  data <- bold::bold_seqspec(taxon = x, sepfasta = FALSE)
  if (length(data) >0 && methods::is(data, "data.frame")) {
    data <- data %>%
      dplyr::na_if("") %>%
      dplyr::filter(markercode == marker) %>% # Remove all sequences for unwanted markers
      dplyr::mutate(domain_name = "Eukaryota") %>%
      dplyr::filter(!is.na(species_name), !is.na(markercode), !is.na(nucleotides)) %>%
      dplyr::filter(!stringr::str_detect(nucleotides, pattern = "I")) #remove inosines

    if (nrow(data) >0) {
      if (output == "standard") {
        data <- data %>%
          dplyr::select(sampleid, species_name, markercode, genbank_accession, nucleotides) %>%
          tidyr::unite("name", c("sampleid", "species_name", "markercode", "genbank_accession"), sep = "|")
      } else if (output == "h") {
        # Hierarchial output
        data <- data %>%
          dplyr::select(sampleid, domain_name, phylum_name, class_name,
                        order_name, family_name, genus_name, species_name,nucleotides
          ) %>%
          tidyr::drop_na() %>%
          tidyr::unite("name", c(
            "sampleid", "domain_name",
            "phylum_name", "class_name",
            "order_name", "family_name",
            "genus_name", "species_name"
          ),
          sep = ";"
          ) %>%
          dplyr::mutate(name = name %>%
                          stringr::str_replace_all(pattern = " ", replacement = "_") %>%
                          trimws(which = "both"))


      } else if (output == "binom") {
        # Binomial output
        data <- data %>%
          dplyr::select(sampleid, species_name, nucleotides) %>%
          tidyr::drop_na() %>%
          tidyr::unite("name", c("sampleid", "species_name"), sep = ";") %>%
          dplyr::mutate(name = name %>%
                          stringr::str_replace_all(pattern = " ", replacement = "_") %>%
                          trimws(which = "both"))


      } else if (output == "bold") {
        # BOLD taxID output
        data <- data %>%
          dplyr::select(sampleid, species_taxID, nucleotides) %>%
          tidyr::drop_na() %>%
          tidyr::unite("name", c("sampleid", "species_taxID"), sep = ";") %>%
          dplyr::mutate(name = name %>%
                          stringr::str_replace_all(pattern = " ", replacement = "_") %>%
                          trimws(which = "both"))

      } else if (output == "gb") {
        # Genbank taxID output
        data <- data %>%
          dplyr::select(sampleid, species_name, nucleotides) %>%
          tidyr::drop_na() %>%
          dplyr::mutate(tax_name = trimws(species_name, which = "both")) %>%
          dplyr::left_join(db, by="tax_name") %>%
          tidyr::unite("name", c("sampleid", "tax_id"), sep = "|")

      } else if (output == "gb-binom") {
        # genbank binomial output
        data <- data %>%
          dplyr::select(sampleid, species_name, nucleotides) %>%
          tidyr::drop_na() %>%
          dplyr::mutate(tax_name = trimws(species_name, which = "both")) %>%
          dplyr::left_join(db, by="tax_name") %>%
          tidyr::unite("name", c("sampleid", "tax_id"), sep = "|") %>%
          tidyr::unite("name", c("name", "tax_name"), sep = ";")
      }

      seqs <- Biostrings::DNAStringSet(data$nucleotides)
      names(seqs) <- data$name
      if (compress == TRUE) {
        Biostrings::writeXStringSet(seqs, out_file, format = "fasta", compress = "gzip", width = 20000)
      } else if (compress == FALSE) {
        Biostrings::writeXStringSet(seqs, out_file, format = "fasta", width = 20000)
      }

      # Done message
      time <- Sys.time() - time
      if (!quiet) (message(paste0("Downloaded ", length(seqs), " ", x, " Sequences from BOLD ", " in ", format(time, digits = 2))))
    }
  }

  # Count number of downloaded sequences
  if(file.exists(out_file)){
    counter <- nrow(Biostrings::fasta.index(out_file))
  } else {
    counter <- 0
  }

  # Count total sequences that should have been downloaded
  if(methods::is(data, "data.frame")){
    total_counter <- nrow(data)
  } else {
    total_counter <- 0
  }

  res <- data.frame(
    taxon = x,
    seqs_total = total_counter,
    seqs_downloaded = counter,
    marker = marker,
    database = "bold",
    time = format(time, digits = 2)
  )
  return(res)
}


#' update_bold
#'
#' @param x A taxon name or vector of taxa to download sequences for
#' @param fasta A fasta file or list of fasta files to check for existing sequence accessions
#' @param marker the barcode marker used as a search term for the database
#' @param quiet Whether progress should be printed to the console.
#' @param output the output format for the taxonomy in fasta headers.
#' Options include "h" for full heirarchial taxonomy (SeqID;Domain;Phylum;Class;Order;Family;Genus;Species),
#' "binom" for just genus species binomials (SeqID;Genus Species),
#' "bold" for BOLD taxonomic ID only (SeqID;BoldTaxID),
#' "gb" for genbank taxonomic ID (SeqID;GBTaxID),
#' "gb-binom" which outputs Genus species binomials, as well as genbank taxonomic ID's, and translates all BOLD taxonomic ID's to genbank taxonomic ID's in the process
#' or "standard" which outputs the default format for each database. For bold this is `sampleid|species name|markercode|genbankid`
#' @param suffix The suffix to add to newly downloaded files. Defaults to 'updates'
#' @param compress Option to compress output fasta files using gzip
#' @param force Option ot overwright files if they already exist
#' @param out_dir Output directory to write fasta files to
#' @param db (Optional) a database file generated using `taxreturn::get_ncbi_taxonomy()` or `taxreturn::get_ott_taxonomy()`
#'
#'
#' @return
#' @import dplyr
#' @import stringr
#' @importFrom bold bold_specimens
#' @importFrom bold bold_seqspec
#' @importFrom tidyr unite
#' @importFrom tidyr drop_na
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings fasta.index
#' @importFrom Biostrings writeXStringSet
#' @importFrom methods is
#'
update_bold <- function(x, fasta, marker = "COI-5P", quiet = FALSE, output = "h", suffix="updates", compress = FALSE, force=FALSE,
                        out_dir = NULL, db=NULL) {

  # function setup
  time <- Sys.time() # get time

  if (!output %in% c("standard", "h", "binom", "gb", "bold", "gb-binom")) {
    stop(paste0(output, " has to be one of: 'standard', 'h','binom','bold', 'gb' or 'gb-binom', see help page for more details"))
  }

  #Check if all inputs exists
  if (!any(file.exists(fasta))) {
    stop("Not all fasta files provided can be found, check the diferectory you have provided")
  }

  #Define directories
  if (is.null(out_dir)) {
    out_dir <- "bold"
    if (!quiet) (message(paste0("No input out_dir given, saving output file to: ", out_dir)))
  }
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }

  # Create output files
  if (compress) {
    out_file <- paste0(normalizePath(out_dir), "/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", marker, "_", suffix, ".fa.gz")
  } else if (!compress) {
    out_file <- paste0(normalizePath(out_dir), "/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", marker, "_", suffix, ".fa")
  }

  #Check if output exists
  if (file.exists(out_file) && force==TRUE) {
    file.remove(out_file)
    cat("", file=out_file)
  } else if (file.exists(out_file) && force==FALSE){
    stop(paste0(out_file, " exists, set force = TRUE to overwrite"))
  }

  # Get NCBI taxonomy database if NCBI format outputs are desired
  if (is.null(db) && output %in% c("gb", "gb-binom")) {
    db <- get_ncbi_taxonomy(include_synonyms = TRUE, force=FALSE)
  }

  # Get accessions from existing fastas
  current <- acc_from_fasta(fasta)

  accs <-  bold::bold_specimens(taxon = x) %>%
    dplyr::pull(sampleid)

  newsearch <- setdiff(accs, current)

  # Bold search
  data <- bold::bold_seqspec(ids = newsearch, sepfasta = FALSE)

  # Find differences

  if (length(data) >0 && !methods::is(data, "logical")) {
    data <- data %>%
      dplyr::na_if("") %>%
      dplyr::filter(stringr::str_detect(marker, markercode)) %>% # Remove all sequences for unwanted markers
      dplyr::mutate(domain_name = "Eukaryota") %>%
      dplyr::filter(!is.na(species_name)) %>%
      dplyr::filter(!stringr::str_detect(nucleotides, pattern = "I")) #remove inosines

    if (nrow(data) >0) {
      if (output == "standard") {
        data <- data %>%
          dplyr::select(sampleid, species_name, markercode, genbank_accession, nucleotides) %>%
          tidyr::unite("name", c("sampleid", "species_name", "markercode", "genbank_accession"), sep = "|")
      } else if (output == "h") {
        # Hierarchial output
        data <- data %>%
          dplyr::select(sampleid, domain_name, phylum_name, class_name,
                        order_name, family_name, genus_name, species_name,nucleotides
          ) %>%
          tidyr::drop_na() %>%
          tidyr::unite("name", c(
            "sampleid", "domain_name",
            "phylum_name", "class_name",
            "order_name", "family_name",
            "genus_name", "species_name"
          ),
          sep = ";"
          ) %>%
          dplyr::mutate(name = name %>%
                          stringr::str_replace_all(pattern = " ", replacement = "_") %>%
                          trimws(which = "both"))


      } else if (output == "binom") {
        # Binomial output
        data <- data %>%
          dplyr::select(sampleid, species_name, nucleotides) %>%
          tidyr::drop_na() %>%
          tidyr::unite("name", c("sampleid", "species_name"), sep = ";") %>%
          dplyr::mutate(name = name %>%
                          stringr::str_replace_all(pattern = " ", replacement = "_") %>%
                          trimws(which = "both"))


      } else if (output == "bold") {
        # BOLD taxID output
        data <- data %>%
          dplyr::select(sampleid, species_taxID, nucleotides) %>%
          tidyr::drop_na() %>%
          tidyr::unite("name", c("sampleid", "species_taxID"), sep = ";") %>%
          dplyr::mutate(name = name %>%
                          stringr::str_replace_all(pattern = " ", replacement = "_") %>%
                          trimws(which = "both"))

      } else if (output == "gb") {
        # Genbank taxID output
        data <- data %>%
          dplyr::select(sampleid, species_name, nucleotides) %>%
          tidyr::drop_na() %>%
          dplyr::mutate(tax_name = trimws(species_name, which = "both")) %>%
          dplyr::left_join(db, by="tax_name") %>%
          tidyr::unite("name", c("sampleid", "tax_id"), sep = "|")

      } else if (output == "gb-binom") {
        # genbank binomial output
        data <- data %>%
          dplyr::select(sampleid, species_name, nucleotides) %>%
          tidyr::drop_na() %>%
          dplyr::mutate(tax_name = trimws(species_name, which = "both")) %>%
          dplyr::left_join(db, by="tax_name") %>%
          tidyr::unite("name", c("sampleid", "tax_id"), sep = "|") %>%
          tidyr::unite("name", c("name", "tax_name"), sep = ";")
      }

      seqs <- Biostrings::DNAStringSet(data$nucleotides)
      names(seqs) <- data$name
      if (compress == TRUE) {
        Biostrings::writeXStringSet(seqs, out_file, format = "fasta", compress = "gzip", width = 20000)
      } else if (compress == FALSE) {
        Biostrings::writeXStringSet(seqs, out_file, format = "fasta", width = 20000)
      }

      # Done message
      time <- Sys.time() - time
      if (!quiet) (message(paste0("Downloaded ", length(seqs), " ", x, " Sequences from BOLD ", " in ", format(time, digits = 2))))
    }
  }
  # Count number of downloaded sequences
  if(file.exists(out_file)){
    counter <- nrow(Biostrings::fasta.index(out_file))
  } else {
    counter <- 0
  }

  # Count total sequences that should have been downloaded
  if(methods::is(data, "data.frame")){
    total_counter <- nrow(data)
  } else {
    total_counter <- 0
  }

  res <- data.frame(
    taxon = x,
    seqs_total = total_counter,
    seqs_downloaded = counter,
    marker = marker,
    database = "bold",
    time = format(time, digits = 2)
  )
  return(res)
}

#' Split bold query
#' This function recursively splits a bold taxonomic query until the amount of records returned is under chunk_size
#'
#' @param x The input taxonomic query
#' @param chunk_size The maximum amount of records to return per query
#' @param quiet Whether progress should be printed to the console.
#'
#' @return
#' @export
#'
#' @import dplyr
#' @importFrom taxize downstream
#' @importFrom bold bold_stats
#' @importFrom methods as
#'
split_bold_query <- function(x, chunk_size=100000, quiet=FALSE){

  rcds <- bold::bold_stats(x, dataType = "overview") %>%
    unlist()
  out <- character()
  if(rcds["total_records"] > chunk_size){
    while(length(x) > 0){
      rcds <- purrr::map(x, ~{
        .x %>%
          bold::bold_stats(dataType = "overview") %>%
          unlist()
      }) %>%
        dplyr::bind_rows() %>%
        dplyr::select("order.count", "family.count", "genus.count", "species.count") %>%
        dplyr::select_if(colSums(.) > nrow(.))

      if(!quiet) {message("Found over ", chunk_size, " (chunk_size value) BOLD records for ", x, ", searching for lower taxonomic ranks to reduce query size")}
      downstream2 <- stringr::str_remove(colnames(rcds[1]), ".count")
      lower_ranks <- taxize::downstream(x, db = "bold", downto = downstream2) %>%
        methods::as("list") %>%
        dplyr::bind_rows() %>%
        dplyr::filter(rank == stringr::str_to_lower(!!downstream2)) %>%
        dplyr::pull(name)

      # Check if any are still over
      rcds2 <- purrr::map(lower_ranks, ~{
        .x %>%
          bold::bold_stats(dataType = "overview") %>%
          unlist()
      }) %>%
        purrr::set_names(lower_ranks) %>%
        dplyr::bind_rows(.id="name")

      # Get successfully resolved
      out <- c(out, rcds2 %>%
                 dplyr::filter(total_records < chunk_size)%>%
                 dplyr::pull(name))

      # Repeat on unresolved
      x <- rcds2 %>%
        dplyr::filter(total_records > chunk_size) %>%
        dplyr::pull(name)
    }

  } else(
    out <- x
  )

  return(out)
}


# wrapper -----------------------------------------------------------------


#' fetch_seqs function
#'
#' @param x A taxon name or vector of taxon names to download sequences for.
#' @param database The database to download from. For NCBI GenBank this currently onlt accepts the arguments 'nuccore' or 'genbank' which is an alias for nuccore.
#' Alternatively sequences can be downloaded from the Barcode of Life Data System (BOLD) using 'bold'
#' @param marker The barcode marker used as a search term for the database. If you are targetting a gene, adding a suffix \[GENE\] will increase the search selectivity.
#' The default for Genbank is 'COI\[GENE\] OR COX1\[GENE\] OR COXI\[GENE\]', while the default for BOLD is 'COI-5P'.
#' If this is set to "mitochondria" and database is 'nuccore', or 'genbank'it will download mitochondrial genomes only.
#' If this is set to "genome" and database is 'nuccore', or 'genbank'it will download complete genome sequences only.
#' @param downstream Instead of search for the query sequence, this provides the option of instead searching for a downstream taxonomic rank.
#' This is useful for big queries where >100k sequences will be downloaded. For example, when x is 'Insecta', and downsteam is Order, this will download all Orders within insecta and thus not overload the query. Default is FALSE.
#' @param output The output format for the taxonomy in fasta headers.
#' Options include "h" for full heirarchial taxonomy (SeqID;Domain;Phylum;Class;Order;Family;Genus;Species),
#' "binom" for just genus species binomials (SeqID;Genus Species),
#' "bold" for BOLD taxonomic ID only (SeqID;BoldTaxID),
#' "gb" for genbank taxonomic ID (SeqID;GBTaxID),
#' "gb-binom" which outputs Genus species binomials, as well as genbank taxonomic ID's, and translates all BOLD taxonomic ID's to genbank taxonomic ID's in the process,
#' or "standard" which outputs the default format for each database. For bold this is `sampleid|species name|markercode|genbankid` while for genbank this is `Accession Sequence definition`
#' @param min_length The maximum length of the query sequence to return. Default 1.
#' @param max_length The maximum length of the query sequence to return.
#' This can be useful for ensuring no off-target sequences are returned. Default 2000.
#' @param subsample (Numeric) return a random subsample of sequences from the search.
#' @param out_dir Output directory to write fasta files to
#' @param compress Option to compress output fasta files using gzip
#' @param force Option to overwrite files if they already exist
#' @param chunk_size Split up the queries made (for genbank), or returned records(for BOLD) into chunks of this size to avoid overloading API servers.
#' if left NULL, the default for genbank searches will be 10,000 for regular queries, 1,000 if marker is "mitochondria", and 1 if marker is "genome"
#' For BOLD queries the default is 100,000 returned records
#' @param multithread Whether multithreading should be used, if TRUE the number of cores will be automatically detected, or provided a numeric vector to manually set the number of cores to use
#' Note, the way this is currently implemented, a seperate worker thread is assigned to each taxon, therefore multithreading will only work
#' if x is a vector, or of downstream is being used.
#' @param quiet Whether progress should be printed to the console.
#' @param progress A logical, for whether or not to print a progress bar when multithread is true. Note, this will slow down processing.
#'
#' @import dplyr
#' @import stringr
#' @import purrr
#' @import future
#' @import furrr
#' @importFrom taxize downstream
#' @importFrom methods as
#'
#' @return
#' @export
#'
fetch_seqs <- function(x, database, marker = NULL, downstream = FALSE,
                       output = "h", min_length = 1, max_length = 2000,
                       subsample=FALSE, chunk_size=NULL, out_dir = NULL, compress = TRUE,
                       force=FALSE, multithread = FALSE, quiet = TRUE, progress=FALSE) {

  if(!database %in% c("nuccore", "genbank", "bold")) {
    stop("database is invalid. See help page for more details")
  }

  if (!output %in% c("standard", "h", "binom", "gb", "bold", "gb-binom")) {
    stop(paste0(output, " has to be one of: 'standard', 'h','binom','bold', 'gb' or 'gb-binom', see help page for more details"))
  }

  #stop if subsample and BOLD is true
  if (is.numeric(subsample ) & database=="bold"){ stop("Subsampling is currently not supported for BOLD")}

  # Get NCBI taxonomy database if NCBI format outputs are desired
  if (database == "bold" && output %in% c("gb", "gb-binom")) {
    db <- get_ncbi_taxonomy(include_synonyms = TRUE, force=FALSE)
  }

  #Define directories
  if (is.null(out_dir)) {
    out_dir <- database
    if (!quiet) (message(paste0("No input out_dir given, saving output file to: ", out_dir)))
  }
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
  out_dir <- normalizePath(out_dir)

  # Evaluate downstream
  if (is.character(downstream)) {
    if (!quiet) cat(paste0("Getting downstream taxa to the level of: ", downstream, "\n"))

    taxlist <- taxize::downstream(x, db = switch(database, bold = "bold", genbank = "ncbi", nuccore = "ncbi"),
                                  downto = downstream) %>%
      methods::as("list") %>%
      dplyr::bind_rows() %>%
      dplyr::filter(rank == stringr::str_to_lower(!!downstream)) %>%
      dplyr::mutate(downloaded = FALSE)

    if (nrow(taxlist) > 0) {
      taxon <- switch(database, bold = taxlist$name, genbank = taxlist$childtaxa_name, nuccore = taxlist$childtaxa_name)
    } else {
      (taxon <- x)
    }
    if (!quiet) cat(paste0(length(taxon), " downstream taxa found\n"))
  } else {
    taxon <- x
  }

  # Setup multithreading - only makes sense if downstream = TRUE
  setup_multithread(multithread = multithread, quiet=quiet)

  # Genbank
  if (database %in% c("genbank", "nuccore")) {
    if (subsample==FALSE) {
      message("Downloading from genbank - No subsampling")
      res <-  furrr::future_map_dfr(
        taxon, search_genbank, database = database, marker = marker,
        output = output, min_length = min_length, max_length = max_length,
        compress = compress, chunk_size=chunk_size, out_dir= out_dir,
        force=force, quiet = quiet, .progress = progress)

    } else if (is.numeric(subsample)){
      message("Downloading from genbank - With subsampling")
      res <-  furrr::future_map_dfr(
        taxon, search_genbank_subsample, database = database, marker = marker,
        out_dir = out_dir, output = output, subsample = subsample,
        min_length = min_length, max_length = max_length, chunk_size=chunk_size,
        force=force, compress = compress, quiet = quiet,  .progress = progress)
    }

  } else if (database == "bold") {
    # Split any querys above chunk_size
    if(is.null(chunk_size)){
      chunk_size <- 100000
    }
    bold_taxon <- furrr::future_map(taxon, split_bold_query, chunk_size=chunk_size, quiet=quiet, .progress = progress) %>%
      unlist()

    if(!quiet) {message("Downloading ", length(bold_taxon)," taxa from BOLD")}
    res <- furrr::future_map(
      bold_taxon, seach_bold, marker = marker, db=db,
      out_dir = out_dir, out_file = NULL, output = output,
      compress = compress, quiet = quiet,  .progress = progress, force=force)
  }

  # Explicitly close multisession workers
  future::plan(future::sequential)

  # Return results summary
  return(res %>%
           dplyr::bind_rows())
}

