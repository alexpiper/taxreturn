#' boldSearch
#' @param x A taxon name or vector of taxa to download sequences for
#' @param marker the barcode marker used as a search term for the database
#' @param quiet (Optional) Print text output
#' @param output the output format for the taxonomy in fasta headers.
#' Options include "h" for full heirarchial taxonomy (SeqID;Domain;Phylum;Class;Order;Family;Genus;Species),
#' "binom" for just genus species binomials (SeqID;Genus Species),
#' "bold" for BOLD taxonomic ID only (SeqID;BoldTaxID),
#' "gb" for genbank taxonomic ID (SeqID;GBTaxID),
#' or "gb-binom" which outputs Genus species binomials, as well as genbank taxonomic ID's, and translates all BOLD taxonomic ID's to genbank taxonomic ID's in the process
#' @param out.file The file to write to, if empty it defaults to the search term
#' @param compress Option to compress output fasta files using gzip
#' @param force Option ot overwright files if they already exist
#' @param out.dir Output directory to write fasta files to
#'
#' @import bold
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @import Biostrings
#'
#' @return
#' @export
#'
#' @examples
boldSearch <- function(x, marker = "COI-5P", quiet = FALSE, output = "h",
                       out.file = NULL, compress = FALSE, force=FALSE,
                       out.dir = NULL, db=NULL) {

  # function setup
  time <- Sys.time() # get time

  if (!output %in% c("h", "binom", "gb", "bold", "gb-binom")) {
    stop(paste0(output, " has to be one of: 'h','binom','bold', 'gb' or 'gb-binom', see help page for more details"))
  }

  # Get NCBI taxonomy database if NCBI format outputs are desired
  if (is.null(db) && output %in% c("gb", "gb-binom")) {
    db <- get_ncbi_lineage(synonyms = TRUE, force=FALSE)
  }

  #Define directories
  if (is.null(out.dir)) {
    out.dir <- "bold"
    if (!quiet) (message(paste0("No input out.dir given, saving output file to: ", out.dir)))
  }
  if (!dir.exists(out.dir)) {
    dir.create(out.dir)
  }

  # Create output files
  if (compress) {
    out.file <- paste0(normalizePath(out.dir), "/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", marker, ".fa.gz")
  } else if (!compress) {
    out.file <- paste0(normalizePath(out.dir), "/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", marker, ".fa")
  }

  #Check if output exists
  if (file.exists(out.file) && force==TRUE) {
    file.remove(out.file)
    cat("", file=out.file)
  } else if (file.exists(out.file) && force==FALSE){
    stop(paste0(out.file, " exists, set force = TRUE to overwrite"))
  }

  # Bold search
  data <- bold::bold_seqspec(taxon = x, sepfasta = FALSE)
  if (length(data) >0 && !class(data) == "logical") {
    data <- data %>%
      dplyr::na_if("") %>%
      dplyr::filter(stringr::str_detect(marker, markercode)) %>% # Remove all sequences for unwanted markers
      dplyr::mutate(domain_name = "Eukaryota") %>%
      dplyr::filter(!is.na(species_name)) %>%
      dplyr::filter(!stringr::str_detect(nucleotides, pattern = "I")) #remove inosines

    if (nrow(data) >0) {
      if (output == "h") {
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
        Biostrings::writeXStringSet(seqs, out.file, format = "fasta", compress = "gzip", width = 20000)
      } else if (compress == FALSE) {
        Biostrings::writeXStringSet(seqs, out.file, format = "fasta", width = 20000)
      }

      # Done message
      time <- Sys.time() - time
      if (!quiet) (message(paste0("Downloaded ", length(seqs), " ", x, " Sequences from BOLD ", " in ", format(time, digits = 2))))
    }
  }
  if(file.exists(out.file)){
    counter <- nrow(Biostrings::fasta.index(out.file))
  } else {
    counter <- 0
  }
  res <- data.frame(
    taxon = x,
    seqs_total = counter,
    seqs_downloaded = counter,
    marker = marker,
    database = "bold",
    time = format(time, digits = 2)
  )
  return(res)
}



# Bold update -------------------------------------------------------------

##  BOLD UPDATE
#' boldUpdate
#'
#' @param x A taxon name or vector of taxa to download sequences for
#' @param fasta A fasta file or list of fasta files to check for existing sequence accessions
#' @param marker the barcode marker used as a search term for the database
#' @param quiet (Optional) Print text output
#' @param output the output format for the taxonomy in fasta headers.
#' Options include "h" for full heirarchial taxonomy (SeqID;Domain;Phylum;Class;Order;Family;Genus;Species),
#' "binom" for just genus species binomials (SeqID;Genus Species),
#' "bold" for BOLD taxonomic ID only (SeqID;BoldTaxID),
#' "gb" for genbank taxonomic ID (SeqID;GBTaxID),
#' or "gb-binom" which outputs Genus species binomials, as well as genbank taxonomic ID's, and translates all BOLD taxonomic ID's to genbank taxonomic ID's in the process
#' @param suffix The suffix to add to newly downloaded files. Defaults to 'updates'
#' @param out.file The file to write to, if empty it defaults to the search term
#' @param compress Option to compress output fasta files using gzip
#' @param force Option ot overwright files if they already exist
#' @param out.dir Output directory to write fasta files to
#' @param db
#'
#'
#' @return
#' @export
#' @import bold
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @import Biostrings
#'
#' @examples
boldUpdate <- function(x, fasta, marker = "COI-5P", quiet = FALSE, output = "h", suffix="updates",
                       out.file = NULL, compress = FALSE, force=FALSE,
                       out.dir = NULL, db=NULL) {

  # function setup
  time <- Sys.time() # get time

  if (!output %in% c("h", "binom", "gb", "bold", "gb-binom")) {
    stop(paste0(output, " has to be one of: 'h','binom','bold', 'gb' or 'gb-binom', see help page for more details"))
  }

  # Get NCBI taxonomy database if NCBI format outputs are desired
  if (is.null(db) && output %in% c("gb", "gb-binom")) {
    db <- get_ncbi_lineage(synonyms = TRUE, force=FALSE)
  }

  #Define directories
  if (is.null(out.dir)) {
    out.dir <- "bold"
    if (!quiet) (message(paste0("No input out.dir given, saving output file to: ", out.dir)))
  }
  if (!dir.exists(out.dir)) {
    dir.create(out.dir)
  }

  # Create output files
  if (compress) {
    out.file <- paste0(normalizePath(out.dir), "/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", marker, "_", suffix, ".fa.gz")
  } else if (!compress) {
    out.file <- paste0(normalizePath(out.dir), "/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", marker, "_", suffix, ".fa")
  }

  #Check if output exists
  if (file.exists(out.file) && force==TRUE) {
    file.remove(out.file)
    cat("", file=out.file)
  } else if (file.exists(out.file) && force==FALSE){
    stop(paste0(out.file, " exists, set force = TRUE to overwrite"))
  }

  # Get accessions from existing fastas
  current <- acc_from_fasta(fasta)

  accs <-  bold::bold_specimens(taxon = x) %>%
    pull(sampleid)

  newsearch <- setdiff(accs, current)

  # Bold search
  data <- bold::bold_seqspec(ids = newsearch, sepfasta = FALSE)

  # Find differences

  if (length(data) >0 && !class(data) == "logical") {
    data <- data %>%
      dplyr::na_if("") %>%
      dplyr::filter(stringr::str_detect(marker, markercode)) %>% # Remove all sequences for unwanted markers
      dplyr::mutate(domain_name = "Eukaryota") %>%
      dplyr::filter(!is.na(species_name)) %>%
      dplyr::filter(!stringr::str_detect(nucleotides, pattern = "I")) #remove inosines

    if (nrow(data) >0) {
      if (output == "h") {
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
        Biostrings::writeXStringSet(seqs, out.file, format = "fasta", compress = "gzip", width = 20000)
      } else if (compress == FALSE) {
        Biostrings::writeXStringSet(seqs, out.file, format = "fasta", width = 20000)
      }

      # Done message
      time <- Sys.time() - time
      if (!quiet) (message(paste0("Downloaded ", length(seqs), " ", x, " Sequences from BOLD ", " in ", format(time, digits = 2))))
    }
  }
  if(file.exists(out.file)){
    counter <- nrow(Biostrings::fasta.index(out.file))
  } else {
    counter <- 0
  }
  res <- data.frame(
    taxon = x,
    seqs_total = counter,
    seqs_downloaded = counter,
    marker = marker,
    database = "bold",
    time = format(time, digits = 2)
  )
  return(res)
}



# Genbank fetching function -----------------------------------------------

#' Genbank search function
#'
#' @param x A taxon name or vector of taxa to download sequences for
#' @param marker The barcode marker used as a search term for the database. If this is set to "mitochondria" it will download full mitochondrial genomes.
#' @param quiet (Optional) Print text output
#' @param output The output format for the taxonomy in fasta headers.
#' Options include "h" for full heirarchial taxonomy (SeqID;Domain;Phylum;Class;Order;Family;Genus;Species),
#' "binom" for just genus species binomials (SeqID;Genus Species),
#' "bold" for BOLD taxonomic ID only (SeqID;BoldTaxID),
#' "gb" for genbank taxonomic ID (SeqID;GBTaxID),
#' or "gb-binom" which outputs Genus species binomials, as well as genbank taxonomic ID's, and translates all BOLD taxonomic ID's to genbank taxonomic ID's in the process
#' @param minlength The minimum length of sequences to download
#' @param maxlength The maximum length of sequences to download
#' @param compress  Option to compress output fasta files using gzip
#' @param force Option ot overwright files if they already exist
#' @param out.dir Output directory to write fasta files to
#'
#'
#' @import rentrez
#' @import biofiles
#' @import Biostrings
#' @import tidyr
#' @import dplyr
#' @import stringr
#' @import purrr
#' @import tibble
#'
#' @return
#' @export
#'
#' Some useful Entrez queries
#'
#' all [filter] 	Retrieves everthing
#' Specified [property] 	Formal binomial and trinomial
#' at or below species level [property]
#' family [rank] 	Rank-based query
#' taxonomy genome [filter] 	Taxa with a direct link to a genome sequence
#' 2009/10/21:2020 [date] 	Date-bounded query
#' mammalia [subtree] 	All taxa within the Mammalia
#' extinct [property] 	Extinct organisms
#' Terminal [property] 	Terminal nodes in the tree
#' loprovencyclife [filter] 	Entries with LinkOut links to the Encyclopedia of Life
#'
#'  # Return res - Taxa - x of y sequences in n chunks in n secs to output
#' @examples
gbSearch <- function(x, database = "nuccore", marker = c("COI", "CO1", "COX1"), quiet = FALSE, output = "h",
                     minlength = 1, maxlength = 2000, subsample=NULL, chunksize=NULL, out.dir = NULL,
                     compress = FALSE, force=FALSE) {

  # function setup
  time <- Sys.time() # get time

  if(!database %in% c("nuccore", "genbank")){
    stop("database is invalid: only nuccore and genbank is currently supported")
  }

  if (!output %in% c("h", "binom", "gb", "bold", "gb-binom")) {
    stop("output has to be one of: 'h','binom','bold', 'gb' or 'gb-binom', see help page for more details")
  }
  if (!tolower(marker) %in% c("mitochondria", "mitochondrion", "genome")) {
    marker <- paste0(paste(marker, collapse="[GENE] OR "),"[GENE]") %>%
      stringr::str_replace_all(stringr::fixed("[GENE][GENE]"), "[GENE]")
  }
  #Define directories
  if (is.null(out.dir)) {
    out.dir <- database
    if (!quiet) (message(paste0("No input out.dir given, saving output file to: ", out.dir)))
  }
  if (!dir.exists(out.dir)) {
    dir.create(out.dir)
  }

  # Create output file
  name <- marker %>%
     stringr::str_replace_all(pattern="\\[GENE]", replacement ="") %>%
     stringr::str_replace_all(pattern="OR ", replacement ="") %>%
     stringr::str_replace_all(pattern=" ", replacement ="_")

  if (compress) {
    out.file <- paste0(normalizePath(out.dir), "/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", name, ".fa.gz")
  } else if (!compress) {
    out.file <- paste0(normalizePath(out.dir), "/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", name, ".fa")
  }

  #Check if output exists
  if (file.exists(out.file) && force==TRUE) {
    file.remove(out.file)
    cat("", file=out.file)
  } else if (file.exists(out.file) && force==FALSE){
    stop(paste0(out.file, " exists, set force = TRUE to overwrite"))
  }

  if(database=="genbank"){
    database="nuccore"
  }

  # Main function
  tryCatch(
    {
      if (!tolower(marker) %in% c("mitochondria", "mitochondrion", "genome")) {
        searchQ <- paste("(", x, "[ORGN])", " AND (", paste(c(marker), collapse = " OR "), ") AND ", minlength, ":", maxlength, "[Sequence Length]", sep = "")
        if(is.null(chunksize)) {chunksize=10000}
      } else if (tolower(marker) %in% c("mitochondria", "mitochondrion")) {
        message(paste0("Input marker is ", marker, ", Downloading full mitochondrial genomes"))
        searchQ <- paste("(", x, "[ORGN])", " AND mitochondrion[filter] AND genome", sep = "")
        if(is.null(chunksize)) {chunksize=1000}
      } else if (tolower(marker) %in% c("genome", "mitochondrion")){
        message(paste0("Input marker is ", marker, ", Downloading full genomes"))
        searchQ <- paste("(", x, "[ORGN])", " AND complete genome[title]", sep = "")
        if(is.null(chunksize)) {chunksize=1}
      }

      search_results <- rentrez::entrez_search(db = database, term = searchQ, retmax = 9999999, use_history = TRUE)

      if (search_results$count > 0 & !is.na(search_results$count)) {
        if (!quiet) {message(paste0(search_results$count, " Sequences to be downloaded for: ", searchQ))}

        l <- 1
        start <- 0

        chunks <- length(search_results$ids) / chunksize
        if (!is.integer(chunks)) {
          chunks <- as.integer(length(search_results$ids) / chunksize) + 1
        }

        for (l in 1:chunks) {

          # Fetch gb flat files
          dl <- rentrez::entrez_fetch(db = database, web_history = search_results$web_history, rettype = "gb", retmax = chunksize, retstart = start)
          gb <- biofiles::gbRecord(rcd = textConnection(dl))

          # Hierarchial output
          if (output == "h") {
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
                                 stringr::str_remove(pattern="\\|.*$")
                               %in% names(seqs)]
          if (compress == TRUE) {
            Biostrings::writeXStringSet(seqs, out.file, format = "fasta", compress = "gzip", width = 20000, append = TRUE)
          } else if (compress == FALSE) {
            Biostrings::writeXStringSet(seqs, out.file, format = "fasta", width = 20000, append = TRUE)
          }

          if (!quiet) (message("Chunk", l, " of ", chunks, " downloaded\r"))
          start <- start + chunksize
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
  if(file.exists(out.file)){
    counter <- nrow(Biostrings::fasta.index(out.file))
  } else {
    counter <- 0
  }
  res <- data.frame(
    taxon = x,
    seqs_total = length(search_results$ids),
    seqs_downloaded = counter,
    marker = marker,
    database = database,
    time = format(time, digits = 2)
  )
  return(res)
}


# Cat_acctax --------------------------------------------------------------

#' Helper function for processing taxids
#'
#' @param x
#'
#' @return
#' @import purrr
#' @import biofiles
#'
#'
#' @examples
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

# Genbank fetching function -----------------------------------------------

#' Genbank subsampling
#'
#' @param x A taxon name or vector of taxa to download sequences for
#' @param marker The barcode marker used as a search term for the database
#' @param quiet
#' @param output The output format for the taxonomy in fasta headers.
#' Options include "h" for full heirarchial taxonomy (SeqID;Domain;Phylum;Class;Order;Family;Genus;Species),
#' "binom" for just genus species binomials (SeqID;Genus Species),
#' "bold" for BOLD taxonomic ID only (SeqID;BoldTaxID),
#' "gb" for genbank taxonomic ID (SeqID;GBTaxID),
#' or "gb-binom" which outputs Genus species binomials, as well as genbank taxonomic ID's, and translates all BOLD taxonomic ID's to genbank taxonomic ID's in the process
#' @param minlength The minimum length of sequences to download
#' @param maxlength The maximum length of sequences to download
#' @param compress  Option to compress output fasta files using gzip
#' @param force Option ot overwright files if they already exist
#' @param out.dir Output directory to write fasta files to
#'
#'
#' @import rentrez
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @import purrr
#' @import biofiles
#' @import Biostrings
#'
#' @return
#' @export
#'
#' @examples
gbSearch_subsample <- function(x, database = "nuccore", marker = c("COI", "CO1", "COX1"),
                               quiet = FALSE, output = "h", minlength = 1, maxlength = 2000,
                               subsample=1000, chunksize=300, compress = FALSE, force=FALSE, out.dir = NULL) {
  # function setup
  time <- Sys.time() # get time

  if(!database %in% c("nuccore", "genbank")){
    stop("database is invalid: only nuccore and genbank is currently supported")
  }

  if (!output %in% c("h", "binom", "gb", "bold", "gb-binom")) {
    stop("output has to be one of: 'h','binom','bold', 'gb' or 'gb-binom', see help page for more details")
  }
  if (!tolower(marker) %in% c("mitochondria", "mitochondrion", "genome")) {
    marker <- paste0(paste(marker, collapse="[GENE] OR "),"[GENE]") %>%
      stringr::str_replace_all(stringr::fixed("[GENE][GENE]"), "[GENE]")
  }
  #Define directories
  if (is.null(out.dir)) {
    out.dir <- database
    if (!quiet) (message(paste0("No input out.dir given, saving output file to: ", out.dir)))
  }
  if (!dir.exists(out.dir)) {
    dir.create(out.dir)
  }

  # Create output file
  name <- marker %>%
    stringr::str_replace_all(pattern="\\[GENE]", replacement ="") %>%
    stringr::str_replace_all(pattern="OR ", replacement ="") %>%
    stringr::str_replace_all(pattern=" ", replacement ="_")

  if (compress) {
    out.file <- paste0(normalizePath(out.dir), "/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", name, ".fa.gz")
  } else if (!compress) {
    out.file <- paste0(normalizePath(out.dir), "/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", name, ".fa")
  }

  #Check if output exists
  if (file.exists(out.file) && force==TRUE) {
    file.remove(out.file)
    cat("", file=out.file)
  } else if (file.exists(out.file) && force==FALSE){
    stop(paste0(out.file, " exists, set force = TRUE to overwrite"))
  }

  if(database=="genbank"){
    database="nuccore"
  }

  # Main function
  tryCatch(
    {
      # Genbank Search
      if (!tolower(marker) %in%c("mitochondria", "mitochondrion", "genome")) {
        searchQ <- paste("(", x, "[ORGN])", " AND (", paste(c(marker), collapse = " OR "), ") AND ", minlength, ":", maxlength, "[Sequence Length]", sep = "")
        chunksize=5000
      } else if (tolower(marker) %in% c("mitochondria", "mitochondrion")) {
        message(paste0("Input marker is ", marker, ", Downloading full mitochondrial genomes"))
        searchQ <- paste("(", x, "[ORGN])", " AND mitochondrion[filter] AND genome", sep = "")
        chunksize=500
      } else if (tolower(marker) %in% c("genome", "mitochondrion")){
        message(paste0("Input marker is ", marker, ", Downloading full genomes"))
        searchQ <- paste("(", x, "[ORGN])", " AND complete genome[title]", sep = "")
        chunksize=1
      }

      search_results <- rentrez::entrez_search(db = database, term = searchQ, retmax = 9999999, use_history = TRUE)

      if (search_results$count > 0 & !is.na(search_results$count)) {

        if (!quiet) {message(paste0(subsample, " of: ", search_results$count, " sequences to be downloaded for: ", searchQ))}

        ids <- sample(search_results$ids, subsample )
        chunks <- split(ids, ceiling(seq_along(ids)/chunksize))

        l <- 1

        for (l in 1:length(chunks)) {

          #Upload ids
          upload <- rentrez::entrez_post(db=database, id=chunks[[l]])

          # Fetch gb flat files
          dl <- rentrez::entrez_fetch(db = database, web_history = upload, rettype = "gb", retmax = 10000)
          gb <- biofiles::gbRecord(rcd = textConnection(dl))

          # Hierarchial output
          if (output == "h") {
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

          # Output FASTA
          seqs <- biofiles::getSequence(gb)
          names(seqs) <- names
          if (compress == TRUE) {
            Biostrings::writeXStringSet(seqs, out.file, format = "fasta", compress = "gzip", width = 20000, append = TRUE)
          } else if (compress == FALSE) {
            Biostrings::writeXStringSet(seqs, out.file, format = "fasta", width = 20000, append = TRUE)
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
  if(file.exists(out.file)){
    counter <- nrow(Biostrings::fasta.index(out.file))
  } else {
    counter <- 0
  }
  res <- data.frame(
    taxon = x,
    seqs_total = length(search_results$ids),
    seqs_downloaded = counter,
    marker = marker,
    database = database,
    time = format(time, digits = 2)
  )
  return(res)
}



# Update genbank database -------------------------------------------------


#' GbUpdate
#'
#' @param x A taxon name or vector of taxa to download sequences for
#' @param fasta A fasta file or list of fasta files to check for existing sequence accessions
#' @param marker The barcode marker used as a search term for the database. If this is set to "mitochondria" it will download full mitochondrial genomes.
#' @param quiet (Optional) Print text output
#' @param output The output format for the taxonomy in fasta headers.
#' Options include "h" for full heirarchial taxonomy (SeqID;Domain;Phylum;Class;Order;Family;Genus;Species),
#' "binom" for just genus species binomials (SeqID;Genus Species),
#' "bold" for BOLD taxonomic ID only (SeqID;BoldTaxID),
#' "gb" for genbank taxonomic ID (SeqID;GBTaxID),
#' or "gb-binom" which outputs Genus species binomials, as well as genbank taxonomic ID's, and translates all BOLD taxonomic ID's to genbank taxonomic ID's in the process
#' @param suffix The suffix to add to newly downloaded files. Defaults to 'updates'
#' @param minlength The minimum length of sequences to download
#' @param maxlength The maximum length of sequences to download
#' @param chunksize The size of the chunked searches to conduct.
#' Warning, chunk sizes over 300 can be too big for the NCBI servers.
#' @param compress  Option to compress output fasta files using gzip
#' @param force Option ot overwright files if they already exist
#' @param out.dir Output directory to write fasta files to
#' @param multithread Whether multithreading should be used
#'
#' @return
#' @export
#' @import stringr
#' @import rentrez
#' @import future
#' @import furrr
#' @import biofiles
#' @import Biostrings
#'
#' @examples
gbUpdate <- function(x, fasta, database = "nuccore", marker = c("COI", "CO1", "COX1"), quiet = FALSE, output = "h", suffix="updates",
                     minlength = 1, maxlength = 2000, chunksize=300, out.dir = NULL,
                     compress = FALSE, force=FALSE, multithread = TRUE){

  # function setup
  time <- Sys.time() # get time

  if(!database %in% c("nuccore", "genbank")){
    stop("database is invalid: only nuccore and genbank is currently supported")
  }

  if (!output %in% c("h", "binom", "gb", "bold", "gb-binom")) {
    stop("output has to be one of: 'h','binom','bold', 'gb' or 'gb-binom', see help page for more details")
  }
  if (!tolower(marker) %in% c("mitochondria", "mitochondrion", "genome")) {
    marker <- paste0(paste(marker, collapse="[GENE] OR "),"[GENE]") %>%
      stringr::str_replace_all(stringr::fixed("[GENE][GENE]"), "[GENE]")
  }
  #Define directories
  if (is.null(out.dir)) {
    out.dir <- database
    if (!quiet) (message(paste0("No input out.dir given, saving output file to: ", out.dir)))
  }
  if (!dir.exists(out.dir)) {
    dir.create(out.dir)
  }

  # Create output file
  name <- marker %>%
    stringr::str_replace_all(pattern="\\[GENE]", replacement ="") %>%
    stringr::str_replace_all(pattern="OR ", replacement ="") %>%
    stringr::str_replace_all(pattern=" ", replacement ="_")

  if (compress) {
    out.file <- paste0(normalizePath(out.dir), "/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", name, "_", suffix, ".fa.gz")
  } else if (!compress) {
    out.file <- paste0(normalizePath(out.dir), "/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", name, "_", suffix, ".fa")
  }

  #Check if output exists
  if (file.exists(out.file) && force==TRUE) {
    file.remove(out.file)
    cat("", file=out.file)
  } else if (file.exists(out.file) && force==FALSE){
    stop(paste0(out.file, " exists, set force = TRUE to overwrite"))
  }

  if(database=="genbank"){
    database="nuccore"
  }

  # Get accessions from existing fastas
  current <- acc_from_fasta(fasta)

  # Genbank Search
  if (!tolower(marker) %in%c("mitochondria", "mitochondrion", "genome")) {
    searchQ <- paste("(", x, "[ORGN])", " AND (", paste(c(marker), collapse = " OR "), ") AND ", minlength, ":", maxlength, "[Sequence Length]", sep = "")
    if(is.null(chunksize)) {chunksize=10000}
  } else if (tolower(marker) %in% c("mitochondria", "mitochondrion")) {
    message(paste0("Input marker is ", marker, ", Downloading full mitochondrial genomes"))
    searchQ <- paste("(", x, "[ORGN])", " AND mitochondrion[filter] AND genome", sep = "")
    if(is.null(chunksize)) {chunksize=1000}
  } else if (tolower(marker) %in% c("genome", "mitochondrion")){
    message(paste0("Input marker is ", marker, ", Downloading full genomes"))
    searchQ <- paste("(", x, "[ORGN])", " AND complete genome[title]", sep = "")
    if(is.null(chunksize)) {chunksize=1}
  }

  search_results <- rentrez::entrez_search(db = database, term = searchQ, retmax = 9999999, use_history = TRUE)

  #Convert search result GIDs to accession
  accs <- gid_to_acc(search_results$ids, db = database, chunksize = chunksize,  multithread=TRUE) %>%
    stringr::str_remove(".[0-9]$")

  # Find any missing accessions
  newsearch <- setdiff(accs, current)

  # Download missing accessions

  if (length(newsearch) > 0) {

    if (!quiet) {message(paste0(length(newsearch), " sequences to be downloaded for: ", searchQ))}

    chunks <- split(newsearch, ceiling(seq_along(newsearch)/chunksize))

    # setup multithreading
    if(multithread){
      future::plan(future::multiprocess)
    } else if(!multithread){
      future::plan(future::sequential)
    }

    #Main function
    furrr::future_map(chunks, function(x){
      upload <- rentrez::entrez_post(db=database, id=x)
      dl <- rentrez::entrez_fetch(db = database, web_history = upload, rettype = "gb", retmax = chunksize)
      gb <- biofiles::gbRecord(rcd = textConnection(dl))

      # Hierarchial output
      if (output == "h") {
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

      # Output FASTA
      seqs <- biofiles::getSequence(gb)
      names(seqs) <- names
      if (compress == TRUE) {
        Biostrings::writeXStringSet(seqs, out.file, format = "fasta", compress = "gzip", width = 20000, append = TRUE)
      } else if (compress == FALSE) {
        Biostrings::writeXStringSet(seqs, out.file, format = "fasta", width = 20000, append = TRUE)
      }
      invisible(NULL)
    })

  }
  #Close all workers
  future::plan(future::sequential)

  #Count how many downloaded
  if(file.exists(out.file)){
    counter <- nrow(Biostrings::fasta.index(out.file))
  } else {
    counter <- 0
  }
  res <- data.frame(
    taxon = x,
    seqs_total = length(search_results$ids),
    seqs_downloaded = counter,
    marker = marker,
    database = database,
    time = format(time, digits = 2)
  )
  return(res)
}




# Fetchseqs wrapper function ----------------------------------------------

#' Fetchseqs wrapper function
#'
#' @param x A taxon name or vector of taxon names to download sequences for.
#' @param database The database to download from. For NCBI GenBank this currently onlt accepts the arguments 'nuccore' or 'genbank' which is an alias for nuccore.
#' Alternatively sequences can be downloaded from the Barcode of Life Data System (BOLD) using 'bold'
#' @param marker The barcode marker used as a search term for the database.
#' The default for Genbank is 'COI OR COI OR COX1 OR COXI', while the default for BOLD is 'COI-5P'.
#' If this is set to "mitochondria" and database is 'nuccore', or 'genbank'it will download mitochondrial genomes only.
#' If this is set to "genome" and database is 'nuccore', or 'genbank'it will download complete genome sequences only.
#' @param downstream Instead of search for the query sequence, this provides the option of instead searching for a downstream taxonomic rank.
#' This is useful for big queries where >100k sequences will be downloaded. For example, when x is 'Insecta', and downsteam is Order, this will download all Orders within insecta and thus not overload the query. Default is FALSE.
#' @param output The output format for the taxonomy in fasta headers.
#' Options include "h" for full heirarchial taxonomy (SeqID;Domain;Phylum;Class;Order;Family;Genus;Species),
#' "binom" for just genus species binomials (SeqID;Genus Species),
#' "bold" for BOLD taxonomic ID only (SeqID;BoldTaxID),
#' "gb" for genbank taxonomic ID (SeqID;GBTaxID),
#' or "gb-binom" which outputs Genus species binomials, as well as genbank taxonomic ID's, and translates all BOLD taxonomic ID's to genbank taxonomic ID's in the process
#' @param minlength The maximum length of the query sequence to return. Default 1.
#' @param maxlength The maximum length of the query sequence to return.
#' This can be useful for ensuring no off-target sequences are returned. Default 2000.
#' @param out.dir Output directory to write fasta files to
#' @param compress Option to compress output fasta files using gzip
#' @param force Option to overwrite files if they already exist
#' @param multithread Whether multithreading should be used.
#' Note, the way this is currently implemented, a seperate worker thread is assigned to each taxon, therefore multithreading will only work
#' if x is a vector, or of downstream is being used.
#' @param quiet (Optional) Print text output
#' @param progress A logical, for whether or not to print a progress bar when multithread is true. Note, this will slow down processing.
#'
#' @import bold
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @import rentrez
#' @import purrr
#' @import future
#' @import furrr
#' @import biofiles
#' @import Biostrings
#' @import taxize
#'
#'
#' @return
#' @export
#'
#' @examples
fetchSeqs <- function(x, database, marker = NULL, downstream = FALSE,
                      output = "h", minlength = 1, maxlength = 2000,
                      subsample=FALSE, chunksize=NULL, out.dir = NULL, compress = TRUE,
                      force=FALSE, multithread = TRUE, quiet = TRUE, progress=FALSE, ...) {

  if(!database %in% c("nuccore", "genbank", "bold")) {
    stop("database is invalid. See help page for more details")
  }

  #stop if subsample and BOLD is true
  if (is.numeric(subsample ) & database=="bold"){ stop("Subsampling is currently not supported for BOLD")}

  # Get NCBI taxonomy database if NCBI format outputs are desired
  if (database == "bold" && output %in% c("gb", "gb-binom")) {
    db <- get_ncbi_lineage(synonyms = TRUE, force=FALSE)
  }

  #Define directories
  if (is.null(out.dir)) {
    out.dir <- database
    if (!quiet) (message(paste0("No input out.dir given, saving output file to: ", out.dir)))
  }
  if (!dir.exists(out.dir)) {
    dir.create(out.dir)
  }
  out.dir <- normalizePath(out.dir)

  #Evaluate downstream
  if (is.character(downstream)) {
    if (!quiet) cat(paste0("Getting downstream taxa to the level of: ", downstream, "\n"))
    taxlist <- taxize::downstream(x, db = "ncbi", downto = downstream) %>%
      as("list") %>%
      dplyr::bind_rows() %>%
      dplyr::filter(rank == stringr::str_to_lower(!!downstream)) %>%
      dplyr::mutate(downloaded = FALSE)

    if (nrow(taxlist) > 0) {
      taxon <- taxlist$childtaxa_name
    } else {
      (taxon <- x)
    }
    if (!quiet) cat(paste0(length(taxon), " downstream taxa found\n"))
  } else {
    taxon <- x
  }

  # setup multithreading - only makes sense if downstream = TRUE
  if(multithread){
    future::plan(future::multiprocess)
  } else if(!multithread){
    future::plan(future::sequential)
  }

  # Genbank
  if (database %in% c("genbank", "nuccore")) {
    if (subsample==FALSE) {
      message("Downloading from genbank - No subsampling")
    res <-  furrr::future_map_dfr(
        taxon, gbSearch, database = database, marker = marker,
        output = output, minlength = minlength, maxlength = maxlength,
        compress = compress, chunksize=chunksize, out.dir= out.dir,
        force=force, quiet = TRUE, .progress = progress, ...=...)

    } else if (is.numeric(subsample)){
      message("Downloading from genbank - With subsampling")
    res <-  furrr::future_map_dfr(
        taxon, gbSearch_subsample, database = database, marker = marker,
        out.dir = out.dir, output = output, subsample = subsample,
        minlength = minlength, maxlength = maxlength, chunksize=chunksize,
        force=force, compress = compress, quiet = TRUE,  .progress = progress, ...=...)
    }

  } else if (database == "bold") {
    checks <- bold::bold_tax_name(taxon)
    bold_taxon <- checks$taxon[which(!is.na(checks$taxon))]
    if (!quiet) {message(paste0(length(bold_taxon), " of ", length(taxon), " taxa found to be valid names on BOLD\n"))}

    message("Downloading from BOLD")
    res <- furrr::future_map(
      bold_taxon, boldSearch, marker = marker, db=db,
      out.dir = out.dir, out.file = NULL, output = output,
      compress = compress, quiet = TRUE,  .progress = progress, force=force, ...=...)
  }

  ## Explicitly close multisession workers
  future::plan(future::sequential)

  # Return results summary
  return(res)
}
