#' boldSearch
#' @param x A taxon name or vector of taxa to download sequences for
#' @param marker the barcode marker used as a search term for the database
#' @param quiet
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
boldSearch <- function(x, marker = NULL, quiet = FALSE, output = "h",
                       out.file = NULL, compress = FALSE, force=FALSE,
                       out.dir = NULL, db=NULL) {

  # function setup
  time <- Sys.time() # get time

  if (!output %in% c("h", "binom", "gb", "bold", "gb-binom")) {
    stop(paste0(output, " has to be one of: 'h','binom','bold', 'gb' or 'gb-binom', see help page for more details"))
  }
  if (is.null(marker)) {
    marker <- "COI-5P"
    if (!quiet) (cat("Using default marker 'COI-5P' \n"))
  }
  if (is.null(out.dir)) {
    out.dir <- "bold"
  }
  if (!file.exists(out.dir)) {
    dir.create(out.dir)
  }
  if (is.null(out.file) & compress == FALSE) {
    out.file <- paste0(getwd(), "/", out.dir, "/", x, "_", marker, ".fa")
  } else if (is.null(out.file) & compress == TRUE) {
    out.file <- paste0(getwd(), "/", out.dir, "/", x, "_", marker, ".fa.gz")
  }
  if (!quiet) (message(paste0("No input file given, saving output file to: ", out.file)))

  if (file.exists(out.file) && force==TRUE) {
    file.remove(out.file)
  } else if (file.exists(out.file) && force==FALSE){
    stop(paste0(out.file, " exists, set force = TRUE to overwright"))
  }
  if (stringr::str_detect(out.file, ".gz")) {
    compress <- TRUE
  }

  # Bold search
  data <- bold::bold_seqspec(taxon = x, sepfasta = FALSE)
  if (length(data) != 0 && !class(data) == "logical") {
    data <- data %>%
      dplyr::na_if("") %>%
      dplyr::filter(grepl(marker, markercode)) %>% # Remove all sequences for unwanted markers
      dplyr::mutate(domain_name = "Eukaryota") %>%
      dplyr::filter(!is.na(species_name))

    if (nrow(data) != 0) {
      # Hierarchial output
      if (output == "h") {
        data <- subset(data, select = c(
          "sampleid", "domain_name",
          "phylum_name", "class_name",
          "order_name", "family_name",
          "genus_name", "species_name",
          "nucleotides"
        )) %>%
          na.omit() %>%
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

        # Binomial output
      } else if (output == "binom") {
        data <- subset(data, select = c("sampleid", "species_name", "nucleotides")) %>%
          na.omit() %>%
          tidyr::unite("name", c("sampleid", "species_name"), sep = ";") %>%
          dplyr::mutate(name = name %>%
                          stringr::str_replace_all(pattern = " ", replacement = "_") %>%
                          trimws(which = "both"))

        # BOLD taxID output
      } else if (output == "bold") {
        data <- subset(data, select = c("sampleid", "species_taxID", "nucleotides")) %>%
          na.omit() %>%
          tidyr::unite("name", c("sampleid", "species_taxID"), sep = ";") %>%
          dplyr::mutate(name = name %>%
                          stringr::str_replace_all(pattern = " ", replacement = "_") %>%
                          trimws(which = "both"))

        # Genbank taxID output
      } else if (output == "gb") {

        data <- subset(data, select = c("sampleid", "species_name", "nucleotides")) %>%
          na.omit() %>%
          dplyr::mutate(tax_name = trimws(species_name, which = "both")) %>%
          dplyr::left_join(db, by="tax_name") %>%
          tidyr::unite("name", c("sampleid", "tax_id"), sep = "|")

        # gb-binom
      } else if (output == "gb-binom") {

        data <- data %>%
          dplyr::select(sampleid, species_name, nucleotides) %>%
          na.omit() %>%
          dplyr::mutate(tax_name = trimws(species_name, which = "both")) %>%
          dplyr::left_join(db, by="tax_name") %>%
          tidyr::unite("name", c("sampleid", "tax_id"), sep = "|") %>%
          tidyr::unite("name", c("name", "tax_name"), sep = ";")
      }

      # Problem -some bold sequences contain an Ionisine
      data <- data %>%
        dplyr::filter(!stringr::str_detect(nucleotides, pattern = "I"))

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
  } else {
    warning(paste0("No species level data for ", x, " on bold\n"))
  }
  invisible(NULL)
}



# Genbank fetching function -----------------------------------------------

#' Genbank search function
#'
#' @param x A taxon name or vector of taxa to download sequences for
#' @param marker The barcode marker used as a search term for the database. If this is set to "mitochondria" it will download full mitochondrial genomes.
#' @param quiet
#' @param output The output format for the taxonomy in fasta headers.
#' Options include "h" for full heirarchial taxonomy (SeqID;Domain;Phylum;Class;Order;Family;Genus;Species),
#' "binom" for just genus species binomials (SeqID;Genus Species),
#' "bold" for BOLD taxonomic ID only (SeqID;BoldTaxID),
#' "gb" for genbank taxonomic ID (SeqID;GBTaxID),
#' or "gb-binom" which outputs Genus species binomials, as well as genbank taxonomic ID's, and translates all BOLD taxonomic ID's to genbank taxonomic ID's in the process
#' @param minlength
#' @param maxlength
#' @param out.file The file to write to, if empty it defaults to the search term
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
#'
#' @examples
gbSearch <- function(x, database = "nuccore", marker = c("COI", "CO1", "COX1"), quiet = FALSE, output = "h",
                     minlength = 1, maxlength = 2000, subsample=NULL,
                     out.file = NULL, compress = FALSE, force=FALSE, out.dir = NULL) {

  # function setup
  time <- Sys.time() # get time

  if(!database %in% c("nuccore")) {
    stop("database is invalid: only nuccore is currently supported")
  }

  if (!output %in% c("h", "binom", "gb", "bold", "gb-binom")) {
    stop("output has to be one of: 'h','binom','bold', 'gb' or 'gb-binom', see help page for more details")
  }
  if (is.null(marker)) {
    marker <- "COI[GENE] OR CO1[GENE] OR COX1[GENE]"
    if (!quiet) (cat("Using default marker 'COI[GENE] OR CO1[GENE] OR COX1[GENE]' \n"))
  } else if(!is.null(marker) && !tolower(marker) %in% c("mitochondria", "mitochondrion", "genome")) {
    marker <- paste0(paste(marker, collapse="[GENE] OR "),"[GENE]")
  }

  if (is.null(out.file)) {
    name <- marker %>%
      stringr::str_replace_all(pattern="\\[GENE]", replacement ="") %>%
      stringr::str_replace_all(pattern="OR ", replacement ="") %>%
      stringr::str_replace_all(pattern=" ", replacement ="_")

    if (!is.null(out.dir) && compress == FALSE) {
      out.file <- paste0(getwd(), "/", out.dir, "/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", name, ".fa")
    } else if (!is.null(out.dir) && compress == TRUE) {
      out.file <- paste0(getwd(), "/", out.dir, "/", stringr::str_replace_all(x,pattern=" ", replacement="_"), "_", name, ".fa.gz")
    } else if (is.null(out.dir) && compress == FALSE) {
      out.file <- paste0(getwd(), "/", "genbank/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", name, ".fa")
    } else if (is.null(out.dir) && compress == TRUE) {
      out.file <- paste0(getwd(), "/", "genbank/", stringr::str_replace_all(x,pattern=" ", replacement="_"), "_", name, ".fa.gz")
    }
    if (!quiet) (message(paste0("No output file given, saving output file to: ", out.file)))
    if (!file.exists("genbank")) {
      dir.create("genbank")
    }
  }
  if (file.exists(out.file) && force==TRUE) {
    file.remove(out.file)
    cat("", file=out.file)
  } else if (file.exists(out.file) && force==FALSE){
    stop(paste0(out.file, " exists, set force = TRUE to overwright"))
  }
  if (stringr::str_detect(out.file, ".gz")) {
    compress <- TRUE
  }


  tryCatch(
    {
      # Genbank Search
      if (!tolower(marker) %in%c("mitochondria", "mitochondrion", "genome")) {
        searchQ <- paste("(", x, "[ORGN])", " AND (", paste(c(marker), collapse = " OR "), ") AND ", minlength, ":", maxlength, "[Sequence Length]", sep = "")
        chunksize=10000
      } else if (tolower(marker) %in% c("mitochondria", "mitochondrion")) {
        message(paste0("Input marker is ", marker, ", Downloading full mitochondrial genomes"))
        searchQ <- paste("(", x, "[ORGN])", " AND mitochondrion[filter] AND genome", sep = "")
        chunksize=1000
      } else if (tolower(marker) %in% c("genome", "mitochondrion")){
        message(paste0("Input marker is ", marker, ", Downloading full genomes"))
        searchQ <- paste("(", x, "[ORGN])", " AND complete genome[title]", sep = "")
        chunksize=1
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

          # Define function for processing taxids
          cat_acctax <- function(x) {
            if(length(biofiles::getAccession(gb)) > 1){
              taxid <- purrr::map(x, biofiles::dbxref, "taxon")
              tax_chr <- purrr::map_chr(taxid, ~{
                as.character(.x[1, 1])})
                attributes(tax_chr) <- NULL
                taxout <- paste0(attributes(taxid)$names, "|", tax_chr)
            } else if(length(biofiles::getAccession(gb)) == 1){
                taxout <- paste0(biofiles::getAccession(x), "|", as.character(biofiles::dbxref(x[1], "taxon")))
              }
            return(taxout)
          }

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
          if(is.na(names(seqs))) {
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
            if (!quiet) {
              counter <- nrow(Biostrings::fasta.index(out.file))
              message(paste0("Downloaded ",counter, " of ", length(seqs), " ", x, " Sequences from Genbank ", " in ", format(time, digits = 2)))
              }
          }
        }
      }
    },
    error = function(e) NULL
  )
  invisible(NULL)
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
#' @param minlength
#' @param maxlength
#' @param out.file The file to write to, if empty it defaults to the search term
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
                               subsample=1000, chunk_size=300, out.file = NULL,
                               compress = FALSE, force=FALSE, out.dir = NULL) {
  # function setup
  time <- Sys.time() # get time

  if(!database %in% c("nuccore")) {
    stop("database is invalid: only nuccore is currently supported")
  }

  if (!output %in% c("h", "binom", "gb", "bold", "gb-binom")) {
    stop("output has to be one of: 'h','binom','bold', 'gb' or 'gb-binom', see help page for more details")
  }
  if (is.null(marker)) {
    marker <- "COI[GENE] OR CO1[GENE] OR COX1[GENE]"
    if (!quiet) (cat("Using default marker 'COI[GENE] OR CO1[GENE] OR COX1[GENE]' \n"))
  } else if(!is.null(marker) && !tolower(marker) %in% c("mitochondria", "mitochondrion", "genome")) {
    marker <- paste0(paste(marker, collapse="[GENE] OR "),"[GENE]")
  }

  if (is.null(out.file)) {
    name <- marker %>%
      stringr::str_replace_all(pattern="\\[GENE]", replacement ="") %>%
      stringr::str_replace_all(pattern="OR ", replacement ="") %>%
      stringr::str_replace_all(pattern=" ", replacement ="_")

    if (!is.null(out.dir) && compress == FALSE) {
      out.file <- paste0(getwd(), "/", out.dir, "/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", name, ".fa")
    } else if (!is.null(out.dir) && compress == TRUE) {
      out.file <- paste0(getwd(), "/", out.dir, "/", stringr::str_replace_all(x,pattern=" ", replacement="_"), "_", name, ".fa.gz")
    } else if (is.null(out.dir) && compress == FALSE) {
      out.file <- paste0(getwd(), "/", "genbank/", stringr::str_replace_all(x, pattern=" ", replacement="_"), "_", name, ".fa")
    } else if (is.null(out.dir) && compress == TRUE) {
      out.file <- paste0(getwd(), "/", "genbank/", stringr::str_replace_all(x,pattern=" ", replacement="_"), "_", name, ".fa.gz")
    }
    if (!quiet) (message(paste0("No input file given, saving output file to: ", out.file)))
    if (!file.exists("genbank")) {
      dir.create("genbank")
    }
  }
  if (file.exists(out.file) && force==TRUE) {
    file.remove(out.file)
    cat("", file=out.file)
  } else if (file.exists(out.file) && force==FALSE){
    stop(paste0(out.file, " exists, set force = TRUE to overwright"))
  }
  if (stringr::str_detect(out.file, ".gz")) {
    compress <- TRUE
  }

  tryCatch(
    {
      # Genbank Search
      if (!tolower(marker) %in%c("mitochondria", "mitochondrion", "genome")) {
        searchQ <- paste("(", x, "[ORGN])", " AND (", paste(c(marker), collapse = " OR "), ") AND ", minlength, ":", maxlength, "[Sequence Length]", sep = "")
        chunksize=10000
      } else if (tolower(marker) %in% c("mitochondria", "mitochondrion")) {
        message(paste0("Input marker is ", marker, ", Downloading full mitochondrial genomes"))
        searchQ <- paste("(", x, "[ORGN])", " AND mitochondrion[filter] AND genome", sep = "")
        chunksize=1000
      } else if (tolower(marker) %in% c("genome", "mitochondrion")){
        message(paste0("Input marker is ", marker, ", Downloading full genomes"))
        searchQ <- paste("(", x, "[ORGN])", " AND complete genome[title]", sep = "")
        chunksize=1
      }

      search_results <- rentrez::entrez_search(db = database, term = searchQ, retmax = 9999999, use_history = TRUE)

      if (search_results$count != 0 & !is.na(search_results$count)) {

        if (!quiet) (message(paste0(subsample, " of: ", search_results$count, " sequences to be downloaded for: ", searchQ)))

        ids <- sample(search_results$ids, subsample )
        chunks <- split(ids, ceiling(seq_along(ids)/chunk_size))

        l <- 1

        for (l in 1:length(chunks)) {

          #Upload ids
          upload <- rentrez::entrez_post(db=database, id=chunks[[l]])

          # Fetch gb flat files
          dl <- rentrez::entrez_fetch(db = database, web_history = upload, rettype = "gb", retmax = 10000)
          gb <- biofiles::gbRecord(rcd = textConnection(dl))

          # Define function for processing taxids
          cat_acctax <- function(x) {
            if(length(x) > 1){
              taxid <- purrr::map(x, biofiles::dbxref, "taxon")
              tax_chr <- purrr::map_chr(taxid, function(y) {
                as.character(y[1, 1])
                attributes(tax_chr) <- NULL
                taxout <- paste0(attributes(taxid)$names, "|", tax_chr)
              })
            } else if(length(x) == 1){
              taxout <- paste0(biofiles::getAccession(x), "|", as.character(biofiles::dbxref(x[1], "taxon")))
            }
            return(taxout)
          }

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
            names <- paste0(biofiles::getAccession(x), ";", biofiles::getOrganism(gb))
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
            if (!quiet) (message(paste0("Downloaded ", length(seqs), " ", x, " Sequences from Genbank ", " in ", format(time, digits = 2))))
          }
        }
      }
    },
    error = function(e) NULL
  )
  return(search_results)
}


# Fetchseqs wrapper function ----------------------------------------------

#' Fetchseqs wrapper function
#'
#' @param x A taxon name or vector of taxa to download sequences for
#' @param database The database to download from. For NCBI GenBank this currently onlt accepts the arguments 'nuccore' or 'genbank' which is an alias for nuccore.
#' Alternatively sequences can be downloaded from the Barcode of Life Data System (BOLD) using 'bold'
#' @param marker The barcode marker used as a search term for the database.
#' The default for Genbank is 'COI OR COI OR COX1 OR COXI', while the default for BOLD is 'COI-5P'.
#' If this is set to "mitochondria" and database is 'nuccore', or 'genbank'it will download mitochondrial genomes only.
#' If this is set to "genome" and database is 'nuccore', or 'genbank'it will download complete genome sequences only.
#' @param downstream Instead of search for the query sequence, this provides the option of instead searching for a downstream taxonomic rank.
#' This is useful for big queries where >100k sequences will be downloaded. For example, when x is 'Insecta', and downsteam is Order, this will download all Orders within insecta and thus not overload the query. Default is FALSE.
#' @param quiet (Optional) Print text output
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
#' @param force Option ot overwright files if they already exist
#' @param cores Number of cores to use
#'
#' @import bold
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @import purrr
#' @import rentrez
#' @import parallel
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
                      quiet = TRUE, output = "h", minlength = 1, maxlength = 2000,
                      subsample=FALSE, out.dir = NULL, compress = TRUE, force=FALSE, cores = 1,...) {

  if(!database %in% c("nuccore", "genbank", "bold")) {
    stop("database is invalid. See help page for more details")
  }
  if (database == "genbank"){
    database <- "nuccore"
  }

  #stop if subsample and BOLD is true
  if (is.numeric(subsample ) & database=="bold"){ stop("Subsampling is currently not supported for BOLD")}

  # Get NCBI Db if NCBI outputs are desired
  if (database == "bold" && output %in% c("gb", "gb-binom")) {
    db <- get_ncbi_lineage(synonyms = TRUE, force=FALSE)
  }

  #Setup parralel
  if (cores == 1) {
    para <- FALSE
    stopclustr <- FALSE
  } else if (cores > 1) {
    #check that cores are available
    navailcores <- parallel::detectCores()
    if(cores > navailcores) stop("Number of cores is more than number available")

    if (!quiet) cat("Multithreading with", cores, "cores\n")
    cores <- parallel::makeCluster(cores, outfile = "out.txt")
    junk <- parallel::clusterEvalQ(cores, sapply(c("bold", "dplyr", "stringr", "purrr", "tidyr", "rentrez", "Biostrings", "biofiles"), require, character.only = TRUE))
    para <- TRUE
    stopclustr <- TRUE
  }

  #Define directories
  if (is.null(out.dir)) {
    out.dir <- database
    if (!quiet) (message(paste0("No input out.dir given, saving output file to: ", out.dir)))
  }
  if (!file.exists(out.dir)) {
    dir.create(out.dir)
  }

  #Evaluate downstream
  if (is.character(downstream)) {
    if (!quiet) cat(paste0("Getting downstream taxa to the level of: ", downstream, "\n"))
    taxlist <- taxize::downstream(x, database = "ncbi", downto = downstream) %>%
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

  #If taxon is a vector
  # Genbank Multithread If possible
  if (database %in% c("nuccore")) {
    if (para == TRUE && subsample==FALSE) {
      message("Multithread downloading from genbank - No subsampling")
      parallel::clusterExport(cl = cores, varlist = c("taxon", "gbSearch", "quiet", "out.dir", "output", "minlength", "maxlength", "compress"), envir = environment())
      parallel::parLapply(cores, taxon, gbSearch, database = database, marker = marker,
                          quiet = quiet, out.file = NULL, output = output,
                          minlength = minlength, maxlength = maxlength,
                          compress = compress, force=force)

    } else if (para == TRUE && is.numeric(subsample)){
      message("Multithread downloading from genbank - With subsampling")
      parallel::clusterExport(cl = cores, varlist = c("taxon", "gbSearch", "quiet", "out.dir", "output", "minlength", "maxlength", "compress"), envir = environment())
      parallel::parLapply(cores, taxon, gbSearch_subsample, database = database, marker = marker,
                          quiet = quiet, out.file = NULL, subsample = subsample,
                          output = output, minlength = minlength,
                          maxlength = maxlength, compress = compress, force=force)
    } else if(para == FALSE && subsample==FALSE){
      message("Sequential downloading from genbank - No subsampling")
      lapply(taxon, gbSearch, database = database, marker = marker,
             quiet = quiet, out.dir = out.dir,
             out.file = NULL, output = output,
             minlength = minlength, maxlength = maxlength,
             compress = compress, force=force)
    } else if(para == FALSE && is.numeric(subsample)){
      message("Sequential downloading from genbank - With subsampling")
      lapply(taxon, gbSearch_subsample, database = database, marker = marker, quiet = quiet, out.dir = out.dir,
             out.file = NULL, output = output, subsample = subsample,
             minlength = minlength, maxlength = maxlength, compress = compress, force=force)
    }

  }
  # BOLD Multithread If possible
  if (database == "bold") {
    if (!quiet) cat("Checking validity of taxon names for BOLD search\n")
    checks <- bold::bold_tax_name(taxon)
    bold_taxon <- checks$taxon[which(!is.na(checks$taxon))]
    if (!quiet) cat(paste0(length(bold_taxon), " of ", length(taxon), " taxa found to be valid names on BOLD\n"))

    # Multithread
    bold_taxon <- if (para == TRUE) {
      message("Multithread downloading from BOLD")
      parallel::parLapply(cores, bold_taxon, boldSearch, marker = marker,
                          quiet = quiet, out.file = NULL, output = output,
                          compress = compress, force=force, db=db)
    } else {
      message("Sequential downloading from BOLD")
      lapply(bold_taxon, boldSearch, marker = marker,
             quiet = quiet, out.dir = out.dir, out.file = NULL,
             output = output, compress = compress, force=force, db=db)
    }
  }

  # Close clusters
  if (para & stopclustr) parallel::stopCluster(cores)
  if (!quiet) message("Done\n")

  if (downstream == TRUE) {
    taxlist$downloaded[which(taxlist$childtaxa_name %in% done)] <- TRUE
    return(taxlist)
  } else {
    return(x)
  }
}

