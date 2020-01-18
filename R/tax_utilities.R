# Get_ranked_lineage ------------------------------------------------------

#' Download NCBI taxdump
#'
#' @param db
#' @param synonyms
#'
#' @return
#' @export
#'
#' @examples
get_ranked_lineage <- function(db = "NCBI", synonyms = TRUE, force=FALSE) {
  if (!identical(db, "NCBI")) {
    stop("Only the NCBI taxonomy database is available in this version\n")
  }
  tmp <- tempdir()
  if (force == TRUE | !file.exists(paste0(tmp, "/rankedlineage.dmp")) | !file.exists(paste0(tmp, "/tmp.tar.gz"))) {
    message("Downloading NCBI taxonomy database")
    fn <- "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
    download.file(fn, destfile = paste0(tmp, "/tmp.tar.gz"))
    message("Extracting data\n")
    test <- untar(tarfile = paste0(tmp, "/tmp.tar.gz"), exdir = tmp)
  if (!identical(test, 0L)) {
    stop(cat(test))
  }
  file.remove(paste0(tmp, "/tmp.tar.gz"))
  }
  message("Building data frame\n")

  # Read data frame
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
  return(lin)
}



# ncbi_taxid --------------------------------------------------------------

#' Get ncbi taxid's for a taxon name
#'
#' @param x
#' @param db
#'
#' @return
#' @export
#'
#' @examples
ncbi_taxid <- function(x, db=NULL) {

  if (is.null(db)) { db <- get_ranked_lineage()}
  out <-  as.data.frame(x) %>%
    magrittr::set_colnames("tax_name") %>%
    dplyr::left_join (db, by = "tax_name") %>%
    dplyr::pull(tax_id)

  return(out)
  }


# Propagate taxonomic assignments to species level ------------------------

#' Propagate taxonomy
#'
#' @param tax
#' @param from
#'
#' @return
#' @export
#'
#' @examples
propagate_tax <- function(tax, from = "Family") {
  col.prefix <- substr(colnames(tax), 1, 1) # Assumes named Kingdom, ...

  # Highest level to propagate from
  if (from == "Phylum") (start <- 2)
  if (from == "Class") (start <- 3)
  if (from == "Order") (start <- 4)
  if (from == "Family") (start <- 5)
  if (from == "Genus") (start <- 6)
  if (from == "Species") (start <- 7)

  # Propagate
  for (col in seq(start, ncol(tax))) {
    prop <- is.na(tax[, col]) & !is.na(tax[, col - 1])
    newtax <- tax[prop, col - 1]
    needs.prefix <- !grepl("^[A-z]__", newtax)
    newtax[needs.prefix] <- paste(col.prefix[col - 1], newtax[needs.prefix], sep = "__")
    tax[prop, col] <- newtax
  }
  tax
}


# summarise fasta ---------------------------------------------------------


#' summarise_fasta
#'
#' @param x The location of a fasta file or gzipped fasta file.
#' @param label optional, Add an extra column with a label
#' @param origin optional, a table with sequence id numbers and their database origins
#'
#' @return
#' @export
#'
#' @examples
summarise_fasta <- function(x, label=NULL, origin=NULL) {
  if(is.null(origin)){
    out <- fasta.index(x) %>%
      mutate(taxid = desc %>%
               str_replace(pattern="(;)(.*?)(?=$)", replacement="")  %>%
               str_replace(pattern="(^)(.*?)(?<=\\|)", replacement="")) %>%
      summarise(nseqs = n(),
                nspecies=n_distinct(taxid),
                mean_length = mean(seqlength),
                q0 = quantile(seqlength, probs=0),
                q25 = quantile(seqlength, probs=0.25),
                q50 = quantile(seqlength, probs=0.5), # Q50 is median
                q75 = quantile(seqlength, probs=0.75),
                q100 = quantile(seqlength, probs=1)
      )

  } else if(is.data.frame(origin) | is_tibble(origin)){
    out <- fasta.index(x) %>%
      mutate(taxid = desc %>%
               str_replace(pattern="(;)(.*?)(?=$)", replacement="")  %>%
               str_replace(pattern="(^)(.*?)(?<=\\|)", replacement="")) %>%
      mutate(seqid = desc %>%
               str_replace(pattern="(\\|)(.*?)(?=$)", replacement=""))  %>%
      left_join(origin, by="seqid") %>%
      group_by(origin) %>%
      summarise(nseqs = n(),
                nspecies=n_distinct(taxid),
                mean_length = mean(seqlength),
                q0 = quantile(seqlength, probs=0),
                q25 = quantile(seqlength, probs=0.25),
                q50 = quantile(seqlength, probs=0.5), # Q50 is median
                q75 = quantile(seqlength, probs=0.75),
                q100 = quantile(seqlength, probs=1)
      )
  }
  if(is.character(label)) {
    out <- out %>%
      mutate(label  = label)
  }
  return(out)
}



# Reformat DADA2 Genus ------------------------------------------------------

# make this just run output heirarchy but removing species rank

#' Reformat DADA2 Genus
#'
#' @param x
#' @param quiet
#' @param ranks
#'
#' @return
#' @export
#'
#' @examples
reformat_dada2_gen <- function(x, db = NULL, quiet = FALSE, ranks = NULL, force=FALSE) {
  if (!is.null(ranks)) {
    ranks <- ranks
  } else if (is.null(ranks)) {
    ranks <- c("kingdom", "phylum", "class", "order", "family", "genus")
    message("No ranks supplied, using default ranks for DADA2 assignTaxonomy: kingdom;phylum;class;order;family;genus")
  }
  x <- reformat_heirarchy(x, db = db, quiet = quiet, ranks = ranks)
  return(x)
}


# Reformat DADA2 Species ----------------------------------------------------

#' Reformat DADA2 Species
#'
#' @param x
#' @param quiet
#'
#' @return
#' @export
#'
#' @examples
reformat_dada2_spp <- function(x, quiet = FALSE) {
  time <- Sys.time() # get time
  # Convert to DNAbin
  if (!is(x, "DNAbin")) {
    x <- ape::as.DNAbin(x)
  }
  seqnames <- stringr::str_split_fixed(names(x), pattern = ";", n = 2) %>%
    as_tibble()
  names(x) <- paste(seqnames$V1, seqnames$V2, sep = " ")
  time <- Sys.time() - time
  if (!quiet) (message(paste0("Reformatted annotations for ", length(x), " sequences in ", format(time, digits = 2))))
  return(x)
}



# Reformat heirarchy ------------------------------------------------------

#' Reformat heirarchy
#'
#' @param x
#' @param quiet
#' @param ranks
#'
#' @return
#' @export
#'
#' @examples
reformat_heirarchy <- function(x, db = NULL, quiet = FALSE, ranks = NULL, sppsep = "_", force=FALSE) {
  time <- Sys.time() # get time
  # Convert to DNAbin
  if (!is(x, "DNAbin")) {
    x <- ape::as.DNAbin(x)
  }
  if (is.null(db)) {
    message("No taxonomy database provided, downloading from NCBI")
    db <- get_ranked_lineage(force=force)
  }

  if (!is.null(ranks)) {
    ranks <- ranks
  } else if (is.null(ranks)) {
    ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
    message("No ranks supplied, using default ranks: kingdom;phylum;class;order;family;genus;species")
  }

  # Split current names
  seqnames <- names(x) %>%
    stringr::str_split_fixed(";", n = 2) %>%
    tibble::as_tibble() %>%
    tidyr::separate(col = V1, into = c("acc", "tax_id"), sep = "\\|") %>%
    dplyr::rename(species = V2) %>%
    dplyr::mutate(tax_id = as.numeric(tax_id))

  # Get lineage from taxid
  lineage <- seqnames %>%
    dplyr::left_join (db %>% dplyr::select(-species), by = "tax_id")  %>%
    tidyr::unite(col = V1, c(acc, tax_id), sep = "|") %>%
    tidyr::unite(col = V2, c(!!ranks), sep = ";")

  names(x) <- paste(lineage$V1, lineage$V2, "", sep = ";")
  time <- Sys.time() - time
  if (!quiet) (message(paste0("Reformatted annotations for ", length(x), " sequences in ", format(time, digits = 2))))
  return(x)
}




# Reformat RDP --------------------------------------------------------------


# Output IDTAXA -----------------------------------------------------------




# taxonomy_to_newick ------------------------------------------------------


## data frame to newick

### recursion function
#traverse <- function(a, i, innerl) {
#  if (i < (ncol(df))) {
#    alevelinner <- as.character(unique(df[which(as.character(df[, i]) == a), i + 1]))
#    desc <- NULL
#    if (length(alevelinner) == 1) {
#      (newickout <- traverse(alevelinner, i + 1, innerl))
#    } else {
#      for (b in alevelinner) desc <- c(desc, traverse(b, i + 1, innerl))
#      il <- NULL
#      if (innerl == TRUE) il <- a
#      (newickout <- paste("(", paste(desc, collapse = ","), ")", il, sep = ""))
#    }
#  }
#  else {
#    (newickout <- a)
#  }
#}
#
## data.frame to newick function
#df2newick <- function(df, innerlabel = FALSE) {
#  alevel <- as.character(unique(df[, 1]))
#  newick <- NULL
#  for (x in alevel) newick <- c(newick, traverse(x, 1, innerlabel))
#  (newick <- paste("(", paste(newick, collapse = ","), ");", sep = ""))
#}
#
#df <- data.frame(x = c("A", "A", "B", "B", "B"), y = c("Ab", "Ac", "Ba", "Ba", "Bd"), z = c("Abb", "Acc", "Bad", "Bae", "Bdd"))
#myNewick <- df2newick(df, TRUE)
#
#library(ape)
#mytree <- read.tree(text = myNewick)
#plot(mytree)


##
