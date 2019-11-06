
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
#' @param out.dir Output directory to write fasta files to
#'
#' @import bold
#' @import tidyverse
#' @import rentrez
#' @import aphid
#' @import insect
#' @import biofiles
#' @import Biostrings
#' @import ape
#' @import stringr
#' @import parallel
#'
#' @return
#' @export
#'
#' @examples
boldSearch <- function(x, marker = NULL, quiet = FALSE, output = "h", out.file = NULL, compress = FALSE, out.dir = NULL) {

  # function setup
  time <- Sys.time() # get time

  # Check if taxizedb is installed
  search_or_sql <- "taxizedb" %in% rownames(installed.packages())
  if (search_or_sql == FALSE) {
    message("taxizedb is not installed, using web queries instead")
  }

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
    out.file <- paste0(getwd(), "/", out.dir, "/", x, "_", marker, "_bold.fa")
  } else if (is.null(out.file) & compress == TRUE) {
    out.file <- paste0(getwd(), "/", out.dir, "/", x, "_", marker, "_bold.fa.gz")
  }
  if (!quiet) (message(paste0("No input file given, saving output file to: ", out.file)))

  if (file.exists(out.file)) {
    file.remove(out.file)
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
          dplyr::mutate(species_name = trimws(species_name, which = "both")) %>%
          dplyr::mutate(gb = taxizedb::name2taxid(data$species_name)) %>%
          #dplyr::mutate(gb = ncbi_taxid(data$species_name)) %>%
          tidyr::unite("name", c("sampleid", "gb"), sep = "|")

        # gb-binom
      } else if (output == "gb-binom") {
        data <- subset(data, select = c("sampleid", "species_name", "nucleotides")) %>%
          na.omit() %>%
          dplyr::mutate(species_name = trimws(species_name, which = "both"))

        ids <- taxizedb::name2taxid(data$species_name, out_type = "summary")
        #ids <- ncbi_taxid(data$species_name)
        # Add exception handling for duplicated taxon names
        if (any(duplicated(ids$name_txt))) {
          dupname <- ids$name_txt[ duplicated(ids$name_txt)]
          dup <- ids[ids$name_txt %in% dupname, ] %>%
            dplyr::mutate(correct = FALSE)
          class <- taxizedb::classification(dup$tax_id)

          for (i in 1:length(class)) {
            dup$correct[i] <- any(stringr::str_detect(class[[i]]$name, pattern = "Insecta"))
          }

          filt <- dup$tax_id[which(dup$correct == FALSE)]
          ids <- ids %>%
            dplyr::filter(!tax_id %in% filt) %>%
            dplyr::rename(species_name = name_txt)

          data <- data %>%
            dplyr::left_join(ids, by = "species_name") %>%
            dplyr::rename(gb = tax_id) %>%
            tidyr::unite("name", c("sampleid", "gb"), sep = "|") %>%
            tidyr::unite("name", c("name", "species_name"), sep = ";")
        } else if (!any(duplicated(ids$name_txt))) {
          data <- data %>%
            dplyr::mutate(gb = taxizedb::name2taxid(data$species_name)) %>%
            #dplyr::mutate(gb = ncbi_taxid(data$species_name)) %>%
            tidyr::unite("name", c("sampleid", "gb"), sep = "|") %>%
            tidyr::unite("name", c("name", "species_name"), sep = ";")
        }
      }

      # Output fASTA
      # iupac <- paste0(paste(names(Biostrings::IUPAC_CODE_MAP),collapse="|"),"|-")

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
#' @param out.dir Output directory to write fasta files to
#'
#'
#' @import bold
#' @import tidyverse
#' @import rentrez
#' @import aphid
#' @import insect
#' @import biofiles
#' @import Biostrings
#' @import ape
#' @import stringr
#' @import parallel
#'
#' @return
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
gbSearch <- function(x, marker = NULL, quiet = FALSE, output = "h", minlength = 1, maxlength = 2000, subsample=NULL, out.file = NULL, compress = FALSE, out.dir = NULL) {

  # function setup
  time <- Sys.time() # get time
  # Check if taxizedb is installed
  search_or_sql <- "taxizedb" %in% rownames(installed.packages())
  if (search_or_sql == FALSE) {
    message("taxizedb is not installed, using web queries instead")
  }

  if (!output %in% c("h", "binom", "gb", "bold", "gb-binom")) {
    stop(paste0(output, " has to be one of: 'h','binom','bold', 'gb' or 'gb-binom', see help page for more details"))
  }
  if (is.null(marker)) {
    marker <- "COI OR COI OR COX1 OR COXI"
    if (!quiet) (cat("Using default marker 'COI OR COI OR COX1 OR COXI' \n"))
  }
  if (is.null(out.file)) {
    if (!is.null(out.dir) && compress == FALSE) {
      out.file <- paste0(getwd(), "/", out.dir, "/", x, "_", marker, "_gb.fa")
    } else if (!is.null(out.dir) && compress == TRUE) {
      out.file <- paste0(getwd(), "/", out.dir, "/", x, "_", marker, "_gb.fa.gz")
    } else if (is.null(out.dir) && compress == FALSE) {
      out.file <- paste0(getwd(), "/", "genbank/", x, "_", marker, "_gb.fa")
    } else if (is.null(out.dir) && compress == TRUE) {
      out.file <- paste0(getwd(), "/", "genbank/", x, "_", marker, "_gb.fa.gz")
    }
    if (!quiet) (message(paste0("No input file given, saving output file to: ", out.file)))
    if (!file.exists("genbank")) {
      dir.create("genbank")
    }
  }
  if (file.exists(out.file)) {
    file.remove(out.file)
  }
  if (stringr::str_detect(out.file, ".gz")) {
    compress <- TRUE
  }

  tryCatch(
    {
      # Genbank Search
      if (!marker %in% c("Mitochondria", "mitochondria","Mito", "mito", "mitochondrion", "Mitochondrion")) {
        searchQ <- paste("(", x, "[ORGN])", " AND (", paste(c(marker), collapse = " OR "), ") AND ", minlength, ":", maxlength, "[Sequence Length]", sep = "")

      } else if (marker %in% c("Mitochondria", "mitochondria","Mito", "mito", "mitochondrion", "Mitochondrion")) {
        message(paste0("Input marker is ", marker, ", Downloading full mitochondrial genomes"))
        searchQ <- paste("(", x, "[ORGN])", " AND mitochondrion[filter] AND genome", sep = "")
        }

      search_results <- rentrez::entrez_search(db = "nuccore", term = searchQ, retmax = 9999999, use_history = TRUE)

      if (search_results$count != 0 & !is.na(search_results$count)) {
        if (!quiet) (message(paste0(search_results$count, " Sequences to be downloaded for: ", searchQ)))

        l <- 1
        start <- 0

        chunks <- length(search_results$ids) / 10000
        if (!is.integer(chunks)) {
          chunks <- as.integer(length(search_results$ids) / 10000) + 1
        }

        for (l in 1:chunks) {

          # Fetch gb flat files
          dl <- rentrez::entrez_fetch(db = "nuccore", web_history = search_results$web_history, rettype = "gb", retmax = 10000, retstart = start)
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
            cat_acctax <- function(x) {
              taxid <- purrr::map(x, biofiles::dbxref, "taxon")
              tax_chr <- purrr::map_chr(taxid, function(y) {
                as.character(y[1, 1])
              })
              attributes(tax_chr) <- NULL
              taxout <- paste0(attributes(taxid)$names, "|", tax_chr)
              return(taxout)
            }
            names <- cat_acctax(gb)
          } else if (output == "binom") {
            names <- paste0(names(biofiles::getSequence(gb)), ";", biofiles::getOrganism(gb))
          } else if (output == "gb-binom") {
            cat_acctax <- function(x) {
              taxid <- purrr::map(x, biofiles::dbxref, "taxon")
              tax_chr <- purrr::map_chr(taxid, function(y) {
                as.character(y[1, 1])
              })
              attributes(tax_chr) <- NULL
              taxout <- paste0(attributes(taxid)$names, "|", tax_chr)
              return(taxout)
            }
            names <- paste0(cat_acctax(gb), ";", biofiles::getOrganism(gb))
          }
          # Output FASTA
          seqs <- biofiles::getSequence(gb) # Something failing here - getting much less than the proper amount

          #Check if names match
          names(seqs) <- names[names %>%
                                 str_replace(pattern="(?<=\\|)(?s)(.*$)", replacement="") %>%
                                 str_replace(pattern="\\|", replacement="")
                               %in% names(seqs)]
          if (compress == TRUE) {
            Biostrings::writeXStringSet(seqs, out.file, format = "fasta", compress = "gzip", width = 20000)
          } else if (compress == FALSE) {
            Biostrings::writeXStringSet(seqs, out.file, format = "fasta", width = 20000)
          }

          if (!quiet) (message("Chunk", l, " of ", chunks, " downloaded\r"))
          start <- start + 10000
          Sys.sleep(2.5)
          if (l >= chunks) {
            time <- Sys.time() - time
            if (!quiet) (message(paste0("Downloaded ", length(seqs), " ", x, " Sequences from Genbank ", " in ", format(time, digits = 2))))
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
#' @param out.dir Output directory to write fasta files to
#'
#'
#' @import bold
#' @import tidyverse
#' @import rentrez
#' @import aphid
#' @import insect
#' @import biofiles
#' @import Biostrings
#' @import ape
#' @import stringr
#' @import parallel
#'
#' @return
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
gbSearch_subsample <- function(x, marker = "COI", quiet = FALSE, output = "h", minlength = 1, maxlength = 2000, subsample=1000, chunk_size=300, out.file = NULL, compress = FALSE, out.dir = NULL) {

  # function setup
  time <- Sys.time() # get time
  # Check if taxizedb is installed
  search_or_sql <- "taxizedb" %in% rownames(installed.packages())
  if (search_or_sql == FALSE) {
    message("taxizedb is not installed, using web queries instead")
  }

  if (!output %in% c("h", "binom", "gb", "bold", "gb-binom")) {
    stop(paste0(output, " has to be one of: 'h','binom','bold', 'gb' or 'gb-binom', see help page for more details"))
  }
  if (is.null(marker)) {
    marker <- "COI OR COI OR COX1 OR COXI"
    if (!quiet) (cat("Using default marker 'COI OR COI OR COX1 OR COXI' \n"))
  }
  if (is.null(out.file)) {
    if (!is.null(out.dir) & compress == FALSE) {
      out.file <- paste0(getwd(), "/", out.dir, "/", x, "_", marker, "_gb.fa")
    } else if (!is.null(out.dir) & compress == TRUE) {
      out.file <- paste0(getwd(), "/", out.dir, "/", x, "_", marker, "_gb.fa.gz")
    } else if (is.null(out.dir) & compress == FALSE) {
      out.file <- paste0(getwd(), "/", "genbank/", x, "_", marker, "_gb.fa")
    } else if (is.null(out.dir) & compress == TRUE) {
      out.file <- paste0(getwd(), "/", "genbank/", x, "_", marker, "_gb.fa.gz")
    }
    if (!quiet) (message(paste0("No input file given, saving output file to: ", out.file)))
    if (!file.exists("genbank")) {
      dir.create("genbank")
    }
  }
  if (file.exists(out.file)) {
    file.remove(out.file)
  }
  if (stringr::str_detect(out.file, ".gz")) {
    compress <- TRUE
  }

  tryCatch(
    {
      # Genbank Search
      if (!marker %in% c("Mitochondria", "mitochondria","Mito", "mito", "mitochondrion", "Mitochondrion")) {
        searchQ <- paste("(", x, "[ORGN])", " AND (", paste(c(marker), collapse = " OR "), ") AND ", minlength, ":", maxlength, "[Sequence Length]", sep = "")

      } else if (marker %in% c("Mitochondria", "mitochondria","Mito", "mito", "mitochondrion", "Mitochondrion")) {
        message(paste0("Input marker is ", marker, ", Downloading full mitochondrial genomes"))
        searchQ <- paste("(", x, "[ORGN])", " AND mitochondrion[filter] AND genome", sep = "")
      }

      if (search_results$count != 0 & !is.na(search_results$count)) {

        if (!quiet) (message(paste0(subsample, " of: ", search_results$count, " sequences to be downloaded for: ", searchQ)))

        ids <- sample(search_results$ids, subsample )
        chunks <- split(ids, ceiling(seq_along(ids)/chunk_size))

        l <- 1

        for (l in 1:length(chunks)) {

          #Upload ids
          upload <- rentrez::entrez_post(db="nuccore", id=chunks[[l]])

          # Fetch gb flat files
          dl <- rentrez::entrez_fetch(db = "nuccore", web_history = upload, rettype = "gb", retmax = 10000)
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
            cat_acctax <- function(x) {
              taxid <- purrr::map(x, biofiles::dbxref, "taxon")
              tax_chr <- purrr::map_chr(taxid, function(y) {
                as.character(y[1, 1])
              })
              attributes(tax_chr) <- NULL
              taxout <- paste0(attributes(taxid)$names, "|", tax_chr)
              return(taxout)
            }
            names <- cat_acctax(gb)
          } else if (output == "binom") {
            names <- paste0(names(biofiles::getSequence(gb)), ";", biofiles::getOrganism(gb))
          } else if (output == "gb-binom") {
            cat_acctax <- function(x) {
              taxid <- purrr::map(x, biofiles::dbxref, "taxon")
              tax_chr <- purrr::map_chr(taxid, function(y) {
                as.character(y[1, 1])
              })
              attributes(tax_chr) <- NULL
              taxout <- paste0(attributes(taxid)$names, "|", tax_chr)
              return(taxout)
            }
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
  invisible(NULL)
}


# Fetchseqs wrapper function ----------------------------------------------

#' Fetchseqs wrapper function
#'
#' @param x A taxon name or vector of taxa to download sequences for
#' @param database
#' @param marker The barcode marker used as a search term for the database.
#' The default for Genbank is 'COI OR COI OR COX1 OR COXI', while the default for BOLD is 'COI-5P'.
#' If this is set to "mitochondria" and database = "genbank" it will download full mitochondrial genomes.
#' @param downstream
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
#' @param cores Number of cores to use
#'


#' @import bold
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @import rentrez
#' @import parallel
#' @import aphid
#' @import insect
#' @import biofiles
#' @import Biostrings
#' @import ape
#' @import taxize
#'
#'
#' @return
#' @export
#'
#' @examples
fetchSeqs <- function(x, database, marker = NULL, downstream = FALSE, quiet = TRUE, output = "h", minlength = 1, maxlength = 2000, subsample=FALSE, out.dir = NULL, compress = TRUE, cores = 1,...) {

  # Check if taxizedb is installed
  search_or_sql <- "taxizedb" %in% rownames(installed.packages())
  if (search_or_sql == FALSE) {
    stop("Error - taxizedb is not installed")
  }

  #stop if subsample and BOLD is true
  if (is.numeric(subsample ) & database=="bold"){ stop("Subsampling is currently not supported for BOLD")}

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
      junk <- parallel::clusterEvalQ(cores, sapply(c("bold", "taxizedb", "tidyverse", "rentrez", "Biostrings", "biofiles"), require, character.only = TRUE))
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
    taxlist <- dplyr::bind_rows(taxize::downstream(x, db = "ncbi", downto = downstream)) %>%
      dplyr::filter(rank == stringr::str_to_lower(!!downstream)) %>%
      dplyr::mutate(downloaded = FALSE)

    if (nrow(taxlist) > 0) {
      taxon <- taxlist$childtaxa_name
    } else {
      (taxon <- x)
    }
    if (!quiet) cat(paste0(length(taxon), " downstream taxa found\n"))
  } else {
    (taxon <- x)
  }

  #If taxon is a vector
  # Genbank Multithread If possible
  if (database == "genbank") {
    taxon <- if (para == TRUE && subsample==FALSE) {
      message("Multithreading with genbank - No subsampling")
      parallel::clusterExport(cl = cores, varlist = c("taxon", "gbSearch", "quiet", "out.dir", "output", "minlength", "maxlength", "compress"), envir = environment())
      parallel::parLapply(cores, taxon, gbSearch, marker = marker,
                          quiet = quiet, out.file = NULL, output = output,
                          minlength = minlength, maxlength = maxlength,
                          compress = compress)

      } else if (para == TRUE && is.numeric(subsample)){
        message("Multithreading with genbank - With subsampling")
        parallel::clusterExport(cl = cores, varlist = c("taxon", "gbSearch", "quiet", "out.dir", "output", "minlength", "maxlength", "compress"), envir = environment())
        parallel::parLapply(cores, taxon, gbSearch_subsample, marker = marker,
                          quiet = quiet, out.file = NULL, subsample = subsample,
                          output = output, minlength = minlength,
                          maxlength = maxlength, compress = compress)
    } else if(para == FALSE && subsample==FALSE){
      message("Sequential processing with genbank - No subsampling")
      lapply(taxon, gbSearch, marker = marker,
             quiet = quiet, out.dir = out.dir,
             out.file = NULL, output = output,
             minlength = minlength, maxlength = maxlength,
             compress = compress)
    } else if(para == FALSE && is.numeric(subsample)){
      message("Sequential processing with genbank - With subsampling")
      lapply(taxon, gbSearch_subsample, marker = marker, quiet = quiet, out.dir = out.dir,
            out.file = NULL, output = output, subsample = subsample,
            minlength = minlength, maxlength = maxlength, compress = compress)
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
      message("Multithreading with BOLD")
      parallel::parLapply(cores, bold_taxon, boldSearch, marker = marker, quiet = quiet, out.file = NULL, output = output, compress = compress)
    } else {
      message("Sequential processing with BOLD")
      lapply(bold_taxon, boldSearch, marker = marker, quiet = quiet, out.dir = out.dir, out.file = NULL, output = output, compress = compress)
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
