# Get_ranked_lineage ------------------------------------------------------

#' Title
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
  if (force == TRUE | !file.exists(paste0(tmp, "/rankedlineage.dmp"))) {
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
