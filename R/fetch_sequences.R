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

#' Fetch sequences from genbank
#'
#' @param x A taxon name or vector of taxa to download sequences for
#' @param database The database to download from. For NCBI GenBank this currently onlt accepts the arguments 'nuccore' or 'genbank' which is an alias for nuccore.
#' @param marker The barcode marker used as a search term for the database.
#' If this is set to "mitochondria" or "mitochondrion" it will download full mitochondrial genomes. If set to "genome" it will download entire genomes only.
#' @param output The output format for the taxonomy in fasta headers.
#' Options include "h" for full heirarchial taxonomy (Accession;Domain;Phylum;Class;Order;Family;Genus;Species),
#' "binom" for just genus species binomials (Accession;Genus_species),
#' "bold" for BOLD taxonomic ID only (Accession;BoldTaxID),
#' "gb" for genbank taxonomic ID (Accession|GBTaxID),
#' "gb-binom" which outputs genbank taxonomic ID's and Genus species binomials, translating BOLD taxonomic ID's to genbank in the process (Accession|GBTaxID;Genus_species)
#' or "standard" which outputs the default format for each database. For bold this is `sampleid|species name|markercode|genbankid`
#' @param min_length The minimum length of sequences to download
#' @param max_length The maximum length of sequences to download
#' @param subsample (Numeric) return a random subsample of sequences from the search.
#' @param chunk_size Split up the query into chunks of this size to avoid overloading API servers. if left NULL, the default will be 300
#' @param db a database file generated using `taxreturn::get_ncbi_taxonomy()`. Generated automatically if NULL.
#' @param multithread Whether multithreading should be used, if TRUE the number of cores will be automatically detected, or provided a numeric vector to manually set the number of cores to use
#' @param quiet Whether progress should be printed to the console.
#' @param progress A logical, for whether or not to print a progress bar when multithread is true. Note, this will slow down processing.
#' @param retry_attempt The number of query attempts in case of query failure due to poor internet connection.
#' @param retry_wait How long to wait between query attempts.
#'
#' @return
#'
#' @examples
fetch_genbank <- function(x, database = "nuccore", marker = c("COI[GENE]", "CO1[GENE]", "COX1[GENE]"), output = "h",
                           min_length = 1, max_length = 2000, subsample=FALSE, chunk_size=100, db=NULL,
                           multithread = FALSE, quiet = FALSE, progress=FALSE, retry_attempt=3, retry_wait=5) {
  if(!database %in% c("nuccore", "genbank")){
    stop("database is invalid: only nuccore and genbank is currently supported")
  }
  if (output == "bold"){
    stop("bold output is only valid for searching the bold database")
  }
  if (!output %in% c("standard", "h", "binom", "gb", "gb-binom")) {
    stop(paste0(output, " has to be one of: 'standard', 'h','binom','gb' or 'gb-binom', see help page for more details"))
  }
  if (!any(tolower(marker) %in% c("mitochondria", "mitochondrion", "genome"))) {
    marker <- paste(marker, collapse=" OR ")
  }

  # Get NCBI taxonomy database if NCBI format outputs are desired
  if (is.null(db) && output %in% c("gb", "gb-binom", "h")) {
    db <- get_ncbi_taxonomy(include_synonyms = TRUE, force=FALSE)
  }

  if(database=="genbank"){
    database="nuccore"
  }
  # Setup multithreading
  setup_multithread(multithread, quiet=quiet)

  # Main function
  if (!tolower(marker) %in% c("mitochondria", "mitochondrion", "genome")) {
    searchQ <- paste("(", x, "[ORGN])", " AND (", paste(c(marker), collapse = " OR "), ") AND ", min_length, ":", max_length, "[Sequence Length]", sep = "")
  } else if (tolower(marker) %in% c("mitochondria", "mitochondrion")) {
    message(paste0("Input marker is ", marker, ", Downloading full mitochondrial genomes"))
    searchQ <- paste("(", x, "[ORGN])", " AND mitochondrion[filter] AND genome", sep = "")
  } else if (tolower(marker) %in% c("genome")){
    message(paste0("Input marker is ", marker, ", Downloading full genomes"))
    searchQ <- paste("(", x, "[ORGN])", " AND complete genome[title]", sep = "")
  }
  if(!quiet){message("Searching genbank with query:", searchQ)}

  search_results <- rentrez::entrez_search(db = database, term = searchQ, retmax = 9999999, use_history = TRUE)
  if (search_results$count > 0 & !is.na(search_results$count)) {

    if (is.numeric(subsample)){
      if (!quiet) {message(paste0(subsample, " of: ", search_results$count, " sequences to be downloaded for: ", searchQ))}
      ids <- sample(search_results$ids, subsample )
    } else{
      if (!quiet) {message(paste0(search_results$count, " Sequences to be downloaded for: ", searchQ))}
      ids <- search_results$ids
    }

    # Split query into chunks
    n <- length(ids)
    r <- rep(1:ceiling(n/chunk_size), each=chunk_size)[1:n]
    id_list <- split(ids, r) %>% magrittr::set_names(NULL)

    # Retrieve each chunk
    seqs <- furrr::future_map(
      id_list, read_genbank_chunk, quiet = FALSE, retry_attempt = retry_attempt, retry_wait = retry_wait, .progress = progress)

    failed <- id_list[!sapply(seqs, class)=="DNAbin"]
    seqs <- seqs[sapply(seqs, class)=="DNAbin"]
    seqs <- concat_DNAbin(seqs)

    # Parse attributes
    seq_spp <- attr(seqs, "species") %>% stringr::str_replace_all("_", " ")
    seq_acc <- attr(seqs, "acc") %>% stringr::str_remove(" .*$")

    if (output == "standard") { # Standard output
      names <- paste0(seq_acc, ";", names(seqs))
    } else if (output == "binom") {
      names <- paste0(seq_acc, ";", seq_spp)
    } else if (output == "h") { # Hierarchical output
      names <- tibble::enframe(seq_spp, name=NULL, value="tax_name") %>%
        dplyr::mutate(sampleid = seq_acc) %>%
        dplyr::left_join(db, by="tax_name") %>%
        tidyr::unite("names", c("sampleid", "superkingdom", "phylum", "class", "order", "family", "genus", "tax_name"), sep = ";")%>%
        dplyr::mutate(names = names %>%
                        stringr::str_replace_all(pattern = ";;", replacement = ";")) %>%
        dplyr::pull(names)
    } else if (output == "gb") { # Genbank taxID output
      names <- tibble::enframe(seq_spp, name=NULL, value="tax_name") %>%
        dplyr::mutate(sampleid = seq_acc) %>%
        dplyr::left_join(db, by="tax_name") %>%
        tidyr::unite("name", c("sampleid", "tax_id"), sep = "|") %>%
        dplyr::pull(name)
    } else if (output == "gb-binom") {
      names <- tibble::enframe(seq_spp, name=NULL, value="tax_name") %>%
        dplyr::mutate(sampleid = seq_acc) %>%
        dplyr::left_join(db, by="tax_name") %>%
        tidyr::unite("name", c("sampleid", "tax_id"), sep = "|") %>%
        tidyr::unite("name", c("name", "tax_name"), sep = ";") %>%
        dplyr::pull(name)
    }
    names(seqs) <- names %>% stringr::str_replace_all(" ", "_")
    out <- seqs
    attr(out, "failed") <- failed
  } else {
    if (!quiet)(message("No sequences available for query: ", searchQ))
    out <- NULL
  }
  # Explicitly close multisession workers
  future::plan(future::sequential)
  return(out)
}

#' Read genbank chunk
#'
#' @param gid a vector of GenBank ID's
#' @param quiet Whether progress should be printed to console.
#' @param retry_attempt The number of query attempts in case of query failure due to poor internet connection.
#' @param retry_wait How long to wait between query attempts
#'
#' @return
#'
#' @examples
read_genbank_chunk <- function(gid, quiet = FALSE, retry_attempt=3, retry_wait=5) {
  n_seqs <- length(gid)
  URL <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=", paste(gid, collapse = ","), "&rettype=gb&retmode=text", sep = "")
  gb <- NULL
  seqs <- NULL
  attempt <- 1
  # Download with error handling
  while((is.null(gb) | (length(seqs) < length(gid))) && attempt <= (retry_attempt+1)) {
    gb <- tryCatch({
      scan(file = URL, what = "", sep = "\n", quiet = TRUE)
    }, error = function(e){
      if (!quiet) {cat(paste("Failed attempt ", attempt,"\n"))}
      Sys.sleep(retry_wait)
      NULL
    })
    if(!is.null(gb)){
      seqs <- parse_gb(gb)
    }
    attempt <- attempt + 1
  }
  if(!length(seqs) == length(gid)){
    if(!quiet)warning("length of returned sequences does not match length of query")
  }
  out <- seqs
  attr(out, "query") <- gid
  return(out)
}

#' Parse genbank flat files
#'
#' @param gb A genbank flat file
#'
#' @return
#' @importFrom stringr str_remove_all
#' @importFrom stringr str_detect
#' @examples
parse_gb <- function(gb){
  # Truncate record at last // to avoid broken records
  gb <- tryCatch({
    gb[1:max(which(grepl("^//", gb)))]
  }, error = function(e){
    warning("Failed parsing")
    print(max(which(grepl("^//", gb))))
    return(NULL)
  })
  if(sum(grepl("ACCESSION", gb)) < 1){
    return(NULL)
  }
  start <- which(grepl("^ORIGIN", gb))
  stop <- which(grepl("^//", gb))

  # Check for malformed start and stops
  good_records <- stringr::str_detect(gb[stop-1], "[0-9] [a-z-]")
  stop <- stop[good_records]

  if (length(start) == length(stop)){
    n_seqs <- length(start)
  } else{
    writeLines(gb, "failed.gb")
    stop("Incorrect length of start and stop, dumped failing records to failed.gb")
  }
  seqs <- vector("character", length = n_seqs)
  for (l in 1:n_seqs){
    seqs[[l]] <- toupper(paste(stringr::str_remove_all(gb[(start[l]+1):(stop[l]-1)], "[^A-Za-z]"), collapse=""))
  }
  # Extract relevant metadata
  seq_acc <- gsub("+ACCESSION +", "", grep("ACCESSION", gb, value = TRUE))
  seq_defs <- gsub("+DEFINITION +", "", grep("DEFINITION", gb, value = TRUE))
  seq_spp <- gsub(" ", "_", gsub(" +ORGANISM +", "", grep(" +ORGANISM +", gb, value = TRUE)))

  names(seqs) <- seq_defs[good_records]
  seqs <- char2DNAbin(seqs)
  attr(seqs, "acc") <- seq_acc[good_records]
  attr(seqs, "species") <- seq_spp[good_records]
  return(seqs)
}

# BOLD functions --------------------------------------------------------

#' Fetch sequences from BOLD
#'
#' @param x A taxon name or vector of taxa to download sequences for
#' @param marker the barcode marker used as a search term for the database
#' @param quiet Whether progress should be printed to the console.
#' @param output the output format for the taxonomy in fasta headers.
#' Options include "h" for full heirarchial taxonomy (Accession;Domain;Phylum;Class;Order;Family;Genus;Species),
#' "binom" for just genus species binomials (Accession;Genus_species),
#' "bold" for BOLD taxonomic ID only (Accession;BoldTaxID),
#' "gb" for genbank taxonomic ID (Accession|GBTaxID),
#' "gb-binom" which outputs genbank taxonomic ID's and Genus species binomials, translating BOLD taxonomic ID's to genbank in the process (Accession|GBTaxID;Genus_species)
#' or "standard" which outputs the default format for each database. For genbank this is the description field, and for bold this is `sampleid|species name|markercode|genbankid`
#' @param db (Optional) a database file generated using `taxreturn::get_ncbi_taxonomy()` or `taxreturn::get_ott_taxonomy()`
#' @param retry_attempt The number of query attempts in case of query failure due to poor internet connection.
#' @param retry_wait How long to wait between query attempts
#'
#' @import dplyr
#' @import stringr
#' @importFrom bold bold_seqspec
#' @importFrom tidyr unite
#' @importFrom tidyr drop_na
#' @importFrom methods is
#' @return
#'
#' @examples
fetch_bold <- function(x, marker = "COI-5P", output = "gb-binom", quiet = FALSE, db=NULL, retry_attempt=3, retry_wait=5) {
  # function setup
  if (!output %in% c("h", "standard", "binom", "gb", "bold", "gb-binom")) {
    stop(paste0(output, " has to be one of: 'standard', 'h','binom','bold', 'gb' or 'gb-binom', see help page for more details"))
  }

  # Get NCBI taxonomy database if NCBI format outputs are desired
  if (is.null(db) && output %in% c("gb", "gb-binom")) {
    db <- get_ncbi_taxonomy(include_synonyms = TRUE, force=FALSE)
  }

  # Bold search
  data <- tryCatch({
    read_bold_chunk(taxon = x, quiet = FALSE, retry_attempt = retry_attempt, retry_wait = retry_wait)
  }, error = function(e){
    writeLines(x, "failed.txt")
    stop("An error during data download, wrote failed taxa to failed.txt")
  })

  bckup_seqs <- data
  if (length(data) >0 & methods::is(data, "data.frame")) {

    data <- tryCatch({
      data %>%
        #dplyr::na_if("") %>%
        dplyr::mutate(domain_name = "Eukaryota") %>%
        dplyr::filter(markercode == marker) %>% # Remove all sequences for unwanted markers
        dplyr::filter(!is.na(species_name), !is.na(markercode), !is.na(nucleotides))
    }, error = function(e){
      saveRDS(bckup_seqs, "bold_data.rds")
      stop("An error during first stage of data parsing, dumped intermediate files to bold_data.rds")
    })

    if (nrow(data) > 0) {

      data <- tryCatch({
        if (output == "standard") {
          data %>%
            dplyr::select(sampleid, species_name, markercode, genbank_accession, nucleotides) %>%
            dplyr::mutate(species_name = species_name %>%
                            stringr::str_replace_all(pattern = " ", replacement = "_") %>%
                            trimws(which = "both")) %>%
            tidyr::unite("name", c("sampleid", "species_name", "markercode", "genbank_accession"), sep = "|")
        } else if (output == "h") {
          # Hierarchial output
          data %>%
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
          data %>%
            dplyr::select(sampleid, species_name, nucleotides) %>%
            tidyr::drop_na() %>%
            tidyr::unite("name", c("sampleid", "species_name"), sep = ";") %>%
            dplyr::mutate(name = name %>%
                            stringr::str_replace_all(pattern = " ", replacement = "_") %>%
                            trimws(which = "both"))
        } else if (output == "bold") {
          # BOLD taxID output
          data %>%
            dplyr::select(sampleid, species_taxID, nucleotides) %>%
            tidyr::drop_na() %>%
            tidyr::unite("name", c("sampleid", "species_taxID"), sep = "|") %>%
            dplyr::mutate(name = name %>%
                            stringr::str_replace_all(pattern = " ", replacement = "_") %>%
                            trimws(which = "both"))

        } else if (output == "gb") {
          # Genbank taxID output
          data %>%
            dplyr::select(sampleid, species_name, nucleotides) %>%
            tidyr::drop_na() %>%
            dplyr::mutate(tax_name = trimws(species_name, which = "both")) %>%
            dplyr::left_join(db, by="tax_name")%>%
            dplyr::mutate(tax_name = tax_name %>%
                            stringr::str_replace_all(pattern = " ", replacement = "_") %>%
                            trimws(which = "both")) %>%
            tidyr::unite("name", c("sampleid", "tax_id"), sep = "|")

        } else if (output == "gb-binom") {
          # genbank binomial output
          data %>%
            dplyr::select(sampleid, species_name, nucleotides) %>%
            tidyr::drop_na() %>%
            dplyr::mutate(tax_name = trimws(species_name, which = "both")) %>%
            dplyr::left_join(db, by="tax_name") %>%
            dplyr::mutate(tax_name = tax_name %>%
                            stringr::str_replace_all(pattern = " ", replacement = "_") %>%
                            trimws(which = "both")) %>%
            tidyr::unite("name", c("sampleid", "tax_id"), sep = "|") %>%
            tidyr::unite("name", c("name", "tax_name"), sep = ";")
        }
      }, error = function(e){
        saveRDS(bckup_seqs, "bold_data.rds")
        stop("An error during second stage of data parsing, dumped intermediate files to bold_data.rds")
      })
      seqs <- char2DNAbin(data$nucleotides)
      names(seqs) <- data$name
      out <- seqs
    } else {
      if (!quiet)(message("No sequences available for query: ", x, " and marker: ", marker))
      out <- NULL
    }
  } else {
    if (!quiet)(message("No sequences available for query: ", x, " and marker: ", marker))
    out <- NULL
  }
  return(out)
}

#' read bold chunk
#'
#' @param taxon A taxon name to download sequences for
#' @param gid a vector of GenBank ID's
#' @param quiet Whether progress should be printed to console.
#' @param retry_attempt The number of query attempts in case of query failure due to poor internet connection.
#' @param retry_wait How long to wait between query attempts
#'
#' @return
#' @export
#'
#' @examples
read_bold_chunk <- function(taxon, quiet = FALSE, retry_attempt=3, retry_wait=5) {
  dl <- NULL
  out <- NULL
  attempt <- 1
  # Download with error handling
  while(is.null(dl) && attempt <= (retry_attempt+1)) {
    dl <- tryCatch({
      cli <- crul::HttpClient$new(url = 'https://v4.boldsystems.org/index.php/API_Public/combined')
      out <- cli$get(query = list(taxon = taxon, format = "tsv"))
      out$raise_for_status()
      if (grepl("html", out$response_headers$`content-type`)) {
        stop(out$parse("UTF-8"))
      }
      out
    }, error = function(e){
      if (!quiet) {cat(paste("Failed attempt ", attempt,"\n"))}
      Sys.sleep(retry_wait)
      NULL
    })
    if(!is.null(dl)){
      tt <- paste0(rawToChar(dl$content, multiple = TRUE), collapse = "")
      if (tt == "") {
        return(NULL)
      }
      Encoding(tt) <- "UTF-8"
      if (grepl("Fatal error", tt)) {
        # set dl to null again to loop through again
        dl <- NULL
      } else{
        out <- readr::read_tsv(tt)
      }
    }
    attempt <- attempt + 1
  }
  return(out)
}


#' Split bold query
#' This function recursively splits a bold taxonomic query until the amount of records returned is under chunk_size
#'
#' @param x The input taxonomic query
#' @param chunk_size The maximum amount of records to return per query
#' @param split_if_under whether to split the query down 1 level even if already below chunksize. This can be useful for making suitable queries for multithread searching.
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
split_bold_query <- function(x, chunk_size=100000, split_if_under = FALSE, quiet=FALSE){

  rcds <- bold::bold_stats(x, dataType = "overview") %>%
    unlist()
  out <- character()
  if(rcds["total_records"] > chunk_size){
    do_split <- TRUE
    if(!quiet) {message("Found over ", chunk_size, " (chunk_size value) BOLD records for ", x, ", searching for lower taxonomic ranks to reduce query size")}
  } else if(split_if_under & (rcds["species.count"] > 1)) {
    do_split <- TRUE
  } else {
    do_split <- FALSE
  }
  if(isTRUE(do_split)){
    while(length(x) > 0){
      rcds <- purrr::map(x, ~{
        .x %>%
          bold::bold_stats(dataType = "overview") %>%
          unlist()
      }) %>%
        dplyr::bind_rows() %>%
        dplyr::select("order.count", "family.count", "genus.count", "species.count") %>%
        dplyr::select_if(colSums(.) > nrow(.))

      downstream2 <- stringr::str_remove(colnames(rcds[1]), ".count")
      lower_ranks <- taxize::downstream(x, db = "bold", downto = downstream2) %>%
        methods::as("list") %>%
        dplyr::bind_rows() %>%
        dplyr::filter(rank == stringr::str_to_lower(!!downstream2)) %>%
        dplyr::pull(name)

      # Check if any are still over (Only if rank is higher than species)
      if (!downstream2 == "species"){
        rcds2 <- purrr::map(lower_ranks, ~{
          .x %>%
            bold::bold_stats(dataType = "overview") %>%
            unlist()
        }) %>%
          purrr::set_names(lower_ranks) %>%
          dplyr::bind_rows(.id="name")
        # Add successfully resolved to out
        out <- c(out, rcds2 %>%
                   dplyr::filter(total_records < chunk_size)%>%
                   dplyr::pull(name))
        # Repeat on unresolved
        x <- rcds2 %>%
          dplyr::filter(total_records > chunk_size) %>%
          dplyr::pull(name)
      } else {
        out <- c(out, lower_ranks)
        x <- character()
      }
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
#' @param output The output format for the taxonomy in fasta headers.
#' Options include "h" for full heirarchial taxonomy (Accession;Domain;Phylum;Class;Order;Family;Genus;Species),
#' "binom" for just genus species binomials (Accession;Genus_species),
#' "bold" for BOLD taxonomic ID only (Accession;BoldTaxID),
#' "gb" for genbank taxonomic ID (Accession|GBTaxID),
#' "gb-binom" which outputs genbank taxonomic ID's and Genus species binomials, translating BOLD taxonomic ID's to genbank in the process (Accession|GBTaxID;Genus_species)
#' or "standard" which outputs the default format for each database. For bold this is `sampleid|species name|markercode|genbankid`
#' @param min_length The maximum length of the query sequence to return. Default 1.
#' @param max_length The maximum length of the query sequence to return.
#' This can be useful for ensuring no off-target sequences are returned. Default 2000.
#' @param subsample (Numeric) return a random subsample of sequences from the search.
#' @param chunk_size Split up the queries made (for genbank), or returned records(for BOLD) into chunks of this size to avoid overloading API servers.
#' if left NULL, the default for genbank searches will be 10,000 for regular queries, 1,000 if marker is "mitochondria", and 1 if marker is "genome"
#' For BOLD queries the default is 100,000 returned records
#' @param db a database file generated using `taxreturn::get_ncbi_taxonomy()`. Generated automatically if NULL.
#' @param multithread Whether multithreading should be used, if TRUE the number of cores will be automatically detected, or provided a numeric vector to manually set the number of cores to use
#' @param quiet Whether progress should be printed to the console.
#' @param progress A logical, for whether or not to print a progress bar when multithread is true. Note, this will slow down processing.
#' @param retry_attempt The number of query attempts in case of query failure due to poor internet connection.
#' @param retry_wait How long to wait between query attempts.
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
#' @examples
fetch_seqs <- function(x, database, marker = NULL, output = "gb-binom",
                        min_length = 1, max_length = 2000, subsample=FALSE,
                        chunk_size=NULL, db=NULL, multithread = FALSE,
                        quiet = FALSE, progress=FALSE, retry_attempt=3, retry_wait=5) {
  # function setup
  time <- Sys.time() # get time

  database <- tolower(database)
  if(!database %in% c("nuccore", "genbank", "bold")) {
    stop("database is invalid. See help page for more details")
  }

  if (!output %in% c("standard", "h", "binom", "gb", "bold", "gb-binom")) {
    stop(paste0(output, " has to be one of: 'standard', 'h','binom','bold', 'gb' or 'gb-binom', see help page for more details"))
  }

  #stop if subsample and BOLD is true
  if (is.numeric(subsample ) & database=="bold"){ stop("Subsampling is currently not supported for BOLD")}

  # Get NCBI taxonomy database if NCBI format outputs are desired
  if (is.null(db) & output %in% c("gb", "gb-binom")) {
    db <- get_ncbi_taxonomy(include_synonyms = TRUE, force=FALSE)
  }

  # Make sure x is unique
  x <- unique(x)

  # Genbank
  if (database %in% c("genbank", "nuccore")) {
    if(!quiet) {message("Downloading from genbank")}
    if(is.null(chunk_size)){
      chunk_size <- 100
    }
    res <- purrr::map(x, fetch_genbank, database = database, marker = marker,
                      output = output, min_length = min_length, max_length = max_length,
                      subsample=subsample, chunk_size=chunk_size, quiet = quiet, multithread=multithread,
                      db=db, retry_attempt = retry_attempt, retry_wait = retry_wait, progress = progress)
    res <- res[!sapply(res, is.null)]
    res <- concat_DNAbin(res)
    res
  } else if (database == "bold") {

    setup_multithread(multithread = multithread, quiet=quiet)

    # Split any very large queries to avoid overloading BOLD query
    if(is.null(chunk_size)){
      chunk_size <- 100000
    }
    # Also If multithreading and only 1 taxon provided, split main query into lower level taxonomy to create something to loop across
    bold_taxon <- furrr::future_map(
      x, split_bold_query, chunk_size=chunk_size, quiet=quiet,
      split_if_under = ifelse((isTRUE(multithread) | (is.numeric(multithread) & multithread > 1)) && length(x) == 1,
                              TRUE, FALSE), .options = furrr::furrr_options(seed = TRUE)) %>%
      unlist()

    if(!quiet) {message("Downloading ", length(bold_taxon)," taxa from BOLD")}
    res <- furrr::future_map(
      bold_taxon, fetch_bold, marker = marker, db=db, output = output,
      quiet = quiet, .progress = progress, .options = furrr::furrr_options(seed = TRUE))
    res <- res[!sapply(res, is.null)]
    res <- concat_DNAbin(res)
    res
  }

  # Explicitly close multisession workers
  future::plan(future::sequential)
  time <- Sys.time() - time
  if (!quiet) (message(paste0("Downloaded ", length(res), " ", x, " Sequences from ",database, " in ", format(time, digits = 2))))
  # Return results summary
  return(res)
}



# Update ------------------------------------------------------------------

#These functions are deprecated and are to be removed

#update_genbank <- function(x, fasta, database = "nuccore", marker = c("COI[GENE]", "CO1[GENE]", "COX1[GENE]"), quiet = FALSE, output = "h", suffix="updates",
#                           min_length = 1, max_length = 2000, chunk_size=300, out_dir = NULL,
#                           compress = FALSE, force=FALSE, multithread = FALSE, progress=FALSE){
#
#  # function setup
#  time <- Sys.time() # get time
#
#  if(!database %in% c("nuccore", "genbank")){
#    stop("database is invalid: only nuccore and genbank is currently supported")
#  }
#
#  #Check if all inputs exists
#  if (!any(file.exists(fasta))) {
#    stop("Not all fasta files provided can be found, check the diferectory you have provided")
#  }
#
#  if (!output %in% c("standard", "h", "binom", "gb", "bold", "gb-binom")) {
#    stop(paste0(output, " has to be one of: 'standard', 'h','binom','bold', 'gb' or 'gb-binom', see help page for more details"))
#  }
#  if (!any(tolower(marker) %in% c("mitochondria", "mitochondrion", "genome"))) {
#    marker <- paste(marker, collapse=" OR ")
#  }
#  #Define directories
#  if (is.null(out_dir)) {
#    out_dir <- database
#    if (!quiet) (message(paste0("No input out_dir given, saving output file to: ", out_dir)))
#  }
#  if (!dir.exists(out_dir)) {
#    dir.create(out_dir)
#  }
#
#  # Create output file
#  name <- marker %>%
#    stringr::str_remove_all(pattern="\\[GENE]") %>%
#    stringr::str_remove_all(pattern="OR ") %>%
#    stringr::str_replace_all(pattern=" ", replacement ="_")
#
#  if (compress) {
#    out_file <- paste0(normalizePath(out_dir), "/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", name, "_", suffix, ".fa.gz")
#  } else if (!compress) {
#    out_file <- paste0(normalizePath(out_dir), "/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", name, "_", suffix, ".fa")
#  }
#
#  #Check if output exists
#  if (file.exists(out_file) && force==TRUE) {
#    file.remove(out_file)
#    cat("", file=out_file)
#  } else if (file.exists(out_file) && force==FALSE){
#    stop(paste0(out_file, " exists, set force = TRUE to overwrite"))
#  }
#
#  if(database=="genbank"){
#    database="nuccore"
#  }
#
#  # Get accessions from existing fastas
#  current <- acc_from_fasta(fasta)
#  if(!quiet){message(length(current), " unique accessions in fastas")}
#
#  # Genbank Search
#  if (!tolower(marker) %in%c("mitochondria", "mitochondrion", "genome")) {
#    searchQ <- paste("(", x, "[ORGN])", " AND (", paste(c(marker), collapse = " OR "), ") AND ", min_length, ":", max_length, "[Sequence Length]", sep = "")
#    if(is.null(chunk_size)) {chunk_size=10000}
#  } else if (tolower(marker) %in% c("mitochondria", "mitochondrion")) {
#    message(paste0("Input marker is ", marker, ", Downloading full mitochondrial genomes"))
#    searchQ <- paste("(", x, "[ORGN])", " AND mitochondrion[filter] AND genome", sep = "")
#    if(is.null(chunk_size)) {chunk_size=1000}
#  } else if (tolower(marker) %in% c("genome", "mitochondrion")){
#    message(paste0("Input marker is ", marker, ", Downloading full genomes"))
#    searchQ <- paste("(", x, "[ORGN])", " AND complete genome[title]", sep = "")
#    if(is.null(chunk_size)) {chunk_size=1}
#  }
#  if(!quiet){message("Searching genbank with query:", searchQ)}
#
#  search_results <- rentrez::entrez_search(db = database, term = searchQ, retmax = 9999999, use_history = TRUE)
#  ids <- search_results$ids
#
#  #Convert search result GIDs to accession
#  if(length(ids) > 0 ){
#    accs <- gid_to_acc(ids, database = database, chunk_size = chunk_size,  multithread=multithread, progress = progress) %>%
#      stringr::str_remove(".[0-9]$")
#  } else{
#    stop("search returned no hits")
#  }
#  # Find any missing accessions
#  newsearch <- setdiff(accs, current)
#
#  # Download missing accessions
#
#  if (length(newsearch) > 0) {
#
#    if (!quiet) {message(paste0(length(newsearch), " sequences to be downloaded for: ", searchQ))}
#
#    chunks <- split(newsearch, ceiling(seq_along(newsearch)/chunk_size))
#
#    # setup multithreading
#    setup_multithread(multithread = multithread, quiet=quiet)
#
#    #Main function
#    furrr::future_map(chunks, function(x){
#      upload <- rentrez::entrez_post(db=database, id=x)
#      dl <- rentrez::entrez_fetch(db = database, web_history = upload, rettype = "gb", retmax = chunk_size)
#      gb <- biofiles::gbRecord(rcd = textConnection(dl))
#
#      # Hierarchial output
#      if (output == "standard") {
#        names <- paste0(biofiles::getAccession(gb), " ", biofiles::getDefinition(gb))
#      } else if (output == "h") {
#        lineage <- biofiles::getTaxonomy(gb) %>%
#          str_split_fixed(pattern = ";", n = Inf) %>%
#          trimws(which = "both") %>%
#          as_tibble() %>%
#          dplyr::mutate(Species = biofiles::getOrganism(gb)) %>%
#          dplyr::mutate(Genus = str_replace(V15, pattern = "[.]", replacement = "")) %>%
#          tidyr::unite("names", c("V1", "V2", "V4", "V6", "V10", "V14", "Genus", "Species"), sep = ";") %>%
#          dplyr::mutate(names = str_replace(names, pattern = " ", replacement = "_"))
#        names <- paste0(names(biofiles::getSequence(gb)), ";", lineage$names)
#
#        # Genbank taxID output
#      } else if (output == "gb") {
#        names <- cat_acctax(gb)
#      } else if (output == "binom") {
#        names <- paste0(biofiles::getAccession(gb), ";", biofiles::getOrganism(gb))
#      } else if (output == "gb-binom") {
#        names <- paste0(cat_acctax(gb), ";", biofiles::getOrganism(gb))
#      }
#      seqs <- biofiles::getSequence(gb) # Something failing here - getting much less than the proper amount
#      if(all(is.na(names(seqs)))) {
#        names(seqs) <- biofiles::getAccession(gb)
#      }
#      #Check if names match
#      names(seqs) <- names[names %>%
#                             stringr::str_remove(pattern=";.*$") %>%
#                             stringr::str_remove(pattern="\\|.*$") %>%
#                             stringr::str_remove(" .*$")
#                           %in% names(seqs)]
#      if (compress == TRUE) {
#        Biostrings::writeXStringSet(seqs, out_file, format = "fasta", compress = "gzip", width = 20000, append = TRUE)
#      } else if (compress == FALSE) {
#        Biostrings::writeXStringSet(seqs, out_file, format = "fasta", width = 20000, append = TRUE)
#      }
#      invisible(NULL)
#    }, .progress = progress)
#
#  }
#  #Close all workers
#  future::plan(future::sequential)
#
#  #Count number of downloaded sequences
#  if(file.exists(out_file)){
#    counter <- nrow(Biostrings::fasta.index(out_file))
#  } else {
#    counter <- 0
#  }
#
#  # Count total sequences that should have been downloaded
#  if(exists("search_results") && length(search_results$ids) > 1 ){
#    total_counter <- length(search_results$ids)
#  } else {
#    total_counter <- 0
#  }
#
#  res <- data.frame(
#    taxon = x,
#    seqs_total = total_counter,
#    seqs_downloaded = counter,
#    marker = marker,
#    database = database,
#    time = format(time, digits = 2)
#  )
#  return(res)
#}

#update_bold <- function(x, fasta, marker = "COI-5P", quiet = FALSE, output = "h", suffix="updates", compress = FALSE, force=FALSE,
#                        out_dir = NULL, db=NULL) {
#
#  # function setup
#  time <- Sys.time() # get time
#
#  if (!output %in% c("standard", "h", "binom", "gb", "bold", "gb-binom")) {
#    stop(paste0(output, " has to be one of: 'standard', 'h','binom','bold', 'gb' or 'gb-binom', see help page for more details"))
#  }
#
#  #Check if all inputs exists
#  if (!any(file.exists(fasta))) {
#    stop("Not all fasta files provided can be found, check the diferectory you have provided")
#  }
#
#  #Define directories
#  if (is.null(out_dir)) {
#    out_dir <- "bold"
#    if (!quiet) (message(paste0("No input out_dir given, saving output file to: ", out_dir)))
#  }
#  if (!dir.exists(out_dir)) {
#    dir.create(out_dir)
#  }
#
#  # Create output files
#  if (compress) {
#    out_file <- paste0(normalizePath(out_dir), "/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", marker, "_", suffix, ".fa.gz")
#  } else if (!compress) {
#    out_file <- paste0(normalizePath(out_dir), "/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", marker, "_", suffix, ".fa")
#  }
#
#  #Check if output exists
#  if (file.exists(out_file) && force==TRUE) {
#    file.remove(out_file)
#    cat("", file=out_file)
#  } else if (file.exists(out_file) && force==FALSE){
#    stop(paste0(out_file, " exists, set force = TRUE to overwrite"))
#  }
#
#  # Get NCBI taxonomy database if NCBI format outputs are desired
#  if (is.null(db) && output %in% c("gb", "gb-binom")) {
#    db <- get_ncbi_taxonomy(include_synonyms = TRUE, force=FALSE)
#  }
#
#  # Get accessions from existing fastas
#  current <- acc_from_fasta(fasta)
#
#  accs <-  bold::bold_specimens(taxon = x) %>%
#    dplyr::pull(sampleid)
#
#  newsearch <- setdiff(accs, current)
#
#  # Bold search
#  data <- bold::bold_seqspec(ids = newsearch, sepfasta = FALSE)
#
#  # Find differences
#
#  if (length(data) >0 && !methods::is(data, "logical")) {
#    data <- data %>%
#      dplyr::na_if("") %>%
#      dplyr::filter(stringr::str_detect(marker, markercode)) %>% # Remove all sequences for unwanted markers
#      dplyr::mutate(domain_name = "Eukaryota") %>%
#      dplyr::filter(!is.na(species_name)) %>%
#      dplyr::filter(!stringr::str_detect(nucleotides, pattern = "I")) #remove inosines
#
#    if (nrow(data) >0) {
#      if (output == "standard") {
#        data <- data %>%
#          dplyr::select(sampleid, species_name, markercode, genbank_accession, nucleotides) %>%
#          tidyr::unite("name", c("sampleid", "species_name", "markercode", "genbank_accession"), sep = "|")
#      } else if (output == "h") {
#        # Hierarchial output
#        data <- data %>%
#          dplyr::select(sampleid, domain_name, phylum_name, class_name,
#                        order_name, family_name, genus_name, species_name,nucleotides
#          ) %>%
#          tidyr::drop_na() %>%
#          tidyr::unite("name", c(
#            "sampleid", "domain_name",
#            "phylum_name", "class_name",
#            "order_name", "family_name",
#            "genus_name", "species_name"
#          ),
#          sep = ";"
#          ) %>%
#          dplyr::mutate(name = name %>%
#                          stringr::str_replace_all(pattern = " ", replacement = "_") %>%
#                          trimws(which = "both"))
#
#
#      } else if (output == "binom") {
#        # Binomial output
#        data <- data %>%
#          dplyr::select(sampleid, species_name, nucleotides) %>%
#          tidyr::drop_na() %>%
#          tidyr::unite("name", c("sampleid", "species_name"), sep = ";") %>%
#          dplyr::mutate(name = name %>%
#                          stringr::str_replace_all(pattern = " ", replacement = "_") %>%
#                          trimws(which = "both"))
#
#
#      } else if (output == "bold") {
#        # BOLD taxID output
#        data <- data %>%
#          dplyr::select(sampleid, species_taxID, nucleotides) %>%
#          tidyr::drop_na() %>%
#          tidyr::unite("name", c("sampleid", "species_taxID"), sep = ";") %>%
#          dplyr::mutate(name = name %>%
#                          stringr::str_replace_all(pattern = " ", replacement = "_") %>%
#                          trimws(which = "both"))
#
#      } else if (output == "gb") {
#        # Genbank taxID output
#        data <- data %>%
#          dplyr::select(sampleid, species_name, nucleotides) %>%
#          tidyr::drop_na() %>%
#          dplyr::mutate(tax_name = trimws(species_name, which = "both")) %>%
#          dplyr::left_join(db, by="tax_name") %>%
#          tidyr::unite("name", c("sampleid", "tax_id"), sep = "|")
#
#      } else if (output == "gb-binom") {
#        # genbank binomial output
#        data <- data %>%
#          dplyr::select(sampleid, species_name, nucleotides) %>%
#          tidyr::drop_na() %>%
#          dplyr::mutate(tax_name = trimws(species_name, which = "both")) %>%
#          dplyr::left_join(db, by="tax_name") %>%
#          tidyr::unite("name", c("sampleid", "tax_id"), sep = "|") %>%
#          tidyr::unite("name", c("name", "tax_name"), sep = ";")
#      }
#
#      seqs <- Biostrings::DNAStringSet(data$nucleotides)
#      names(seqs) <- data$name
#      if (compress == TRUE) {
#        Biostrings::writeXStringSet(seqs, out_file, format = "fasta", compress = "gzip", width = 20000)
#      } else if (compress == FALSE) {
#        Biostrings::writeXStringSet(seqs, out_file, format = "fasta", width = 20000)
#      }
#
#      # Done message
#      time <- Sys.time() - time
#      if (!quiet) (message(paste0("Downloaded ", length(seqs), " ", x, " Sequences from BOLD ", " in ", format(time, digits = 2))))
#    }
#  }
#  # Count number of downloaded sequences
#  if(file.exists(out_file)){
#    counter <- nrow(Biostrings::fasta.index(out_file))
#  } else {
#    counter <- 0
#  }
#
#  # Count total sequences that should have been downloaded
#  if(methods::is(data, "data.frame")){
#    total_counter <- nrow(data)
#  } else {
#    total_counter <- 0
#  }
#
#  res <- data.frame(
#    taxon = x,
#    seqs_total = total_counter,
#    seqs_downloaded = counter,
#    marker = marker,
#    database = "bold",
#    time = format(time, digits = 2)
#  )
#  return(res)
#}
