# DNAbin2DNAStringSet -----------------------------------------------------

#' Convert DNABin to DNAStringSet
#'
#' @param x a DNABin object
#' @param remove_gaps Whether gaps should be removed
#'
#' @return
#' @export
#' @import ape
#' @import Biostrings
#' @import purrr
#'
#' @examples
DNAbin2DNAstringset <- function (x, remove_gaps = FALSE) {
  if(!is(x, "DNAbin")){
    stop("Input must be a DNAbin")
  }
  x %>%
    as.list() %>%
    as.character %>%
    purrr::map(function(y){
      y[!y %in% tolower(Biostrings::DNA_ALPHABET)] <- "N"
      if(isTRUE(remove_gaps)){
        y[y=="-"] <- ""
      }
      paste0(y, collapse="")
    })%>%
    unlist %>%
    Biostrings::DNAStringSet()
}


# Accession to hexadecimal coding ----------------------------------------------


#' Convert accession number to hexadecimal coding
#' @description Hexadecimal coding is useful for ensuring that accession numbers dont contain characters that are considered illegal for newick strings
#'
#' @param x A DNAbin or DNAStringSet with names formatted Accession|taxid;taxonomy.
#' Or a character vector with names formatted Accession|taxid;taxonomy
#' Or a character vector of accessions
#'
#' @return
#' @export
#' @import sodium
#' @import ape
#' @import Biostrings
#' @import stringr
#'
#' @examples
acc2hex <- function(x){
  #Define hex function
  .hex <- function(y){
    bin2hex(sodium::charToRaw(y))
  }
  #Check input format
  if (is(x, "DNAbin") | is(x, "DNAStringSet")) {
    message("Input is ", class(x),", assuming header is in 'Accession|taxid;taxonomy' format")
    acc <- names(x) %>%
      stringr::str_remove("\\|.*$") %>%
      purrr::map_chr(.hex)
    names(x) <- names(x) %>%
      stringr::str_remove("^.*(?=\\|)") %>%
      paste0(acc, .)
  }else  if (is(x, "character") && (stringr::str_detect(x, "\\|") & stringr::str_detect(x, ";"))) {
    message("Detected | and ; delimiters, assuming 'Accession|taxid;taxonomy' format")
    acc <- x %>%
      stringr::str_remove("\\|.*$") %>%
      purrr::map_chr(.hex)
    x <- x %>%
      stringr::str_remove("^.*(?=\\|)") %>%
      paste0(acc, .)
  } else  if (is(x, "character") && !(stringr::str_detect(x, "\\|") && stringr::str_detect(x, ";"))) {
    message("Did not detect | and ; delimiters, assuming a vector of accessions")
    x <- x %>%
      purrr::map_chr(.hex)
  }else (stop("x must be DNA bin, DNAStringSet or vector of accesssions"))
  return(x)
}


# Hexadecimal coding to accession -----------------------------------------

#' Title
#' @description Hexadecimal coding is useful for ensuring that accession numbers dont contain characters that are considered illegal for newick strings
#' @param x A DNAbin or DNAStringSet with names formatted Accession|taxid;taxonomy.
#' Or a character vector with names formatted Accession|taxid;taxonomy
#' Or a character vector of accessions
#'
#' @return
#' @export
#' @import sodium
#' @import ape
#' @import Biostrings
#' @import stringr
#'
#' @examples
hex2acc <- function(x){
  #Define hex function
  .unhex <- function(y){
    rawToChar(sodium::hex2bin(y))
  }
  #Check input format
  if (is(x, "DNAbin") | is(x, "DNAStringSet")) {
    message("Input is ", class(x),", assuming header is in 'Accession|taxid;taxonomy' format")
    acc <- names(x) %>%
      stringr::str_remove("\\|.*$") %>%
      purrr::map_chr(.unhex)
    names(x) <- names(x) %>%
      str_remove("^.*(?=\\|)") %>%
      paste0(acc, .)
  }else  if (is(x, "character") && (stringr::str_detect(x, "\\|") & stringr::str_detect(x, ";"))) {
    message("Detected | and ; delimiters, assuming 'Accession|taxid;taxonomy' format")
    acc <- x %>%
      stringr::str_remove("\\|.*$") %>%
      purrr::map_chr(.unhex)
    x <- x %>%
      stringr::str_remove("^.*(?=\\|)") %>%
      paste0(acc, .)
  } else  if (is(x, "character") && !(stringr::str_detect(x, "\\|") && stringr::str_detect(x, ";"))) {
    message("Did not detect | and ; delimiters, assuming a vector of accessions")
    x <- x %>%
      purrr::map_chr(.unhex)
  }else (stop("x must be DNA bin, DNAStringSet or vector of accesssions"))
  return(x)
}
