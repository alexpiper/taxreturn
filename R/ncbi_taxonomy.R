#' Download NCBI taxdump
#'
#' @param db
#' @param synonyms
#'
#' @return
#' @export
#'
#' @examples
get_ncbi_lineage <- function(db = "NCBI", synonyms = TRUE, force=FALSE) {
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
  attr(lin,'type') <- 'ncbi'
  return(lin)
}


#' (Deprecated) Get ranked Lineage
#'
#' @param db
#' @param synonyms
#' @param force
#'
#' @return
#' @export
#'
#' @examples
get_ranked_lineage <- function(db = "NCBI", synonyms = TRUE, force=FALSE) {
  .Deprecated("get_ncbi_lineage") #include a package argument, too
  get_ncbi_lineage(db = "NCBI", synonyms = TRUE, force = FALSE)
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

  if (is.null(db)) { db <- get_ncbi_lineage()}
  out <-  as.data.frame(x) %>%
    magrittr::set_colnames("tax_name") %>%
    dplyr::left_join (db, by = "tax_name") %>%
    dplyr::pull(tax_id)

  return(out)
}
