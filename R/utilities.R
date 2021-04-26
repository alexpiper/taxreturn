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

# Function to generate random sequences
#' Generate random sequences
#'
#' @param n The number of sequences to generate
#' @param length The length of the generated sequences
#' @param alphabet The DNA alphabet to draw from, default is A,G,C,T
#'
#' @return
#' @export
#'
#' @examples
random_seq <- function(n, length, alphabet = c("A","G","T","C")){
  out <- seq(1, n, 1) %>%
    purrr::map2(length, function(x,y ){
      paste(sample(alphabet, y, replace = T), collapse="")
    }) %>%
    insect::char2dna()
  names(out) <- make.unique(rep("Seq", n), sep="_")
  return(out)
}

# Accession to hexadecimal coding ----------------------------------------------


#' Convert accession number to hexadecimal coding
#' @description Hexadecimal coding is useful for ensuring that accession numbers dont contain characters that are considered illegal for newick strings
#'
#' @param x A DNAbin or DNAStringSet with names formatted Accession|taxid;taxonomy.
#' Or a character vector with names formatted Accession|taxid;taxonomy
#' Or a character vector of accessions
#' @param force override checks if string is already hexadecimal
#'
#' @return
#' @export
#' @import ape
#' @import Biostrings
#' @import stringr
#'
#' @examples
acc2hex <- function(x, force=FALSE){
  #Check input format
  if (is(x, "DNAbin") | is(x, "DNAStringSet")) {
    message("Input is ", class(x),", assuming header is in 'Accession|taxid;taxonomy' format")
    #Check if already hexed
    if(.ishex(names(x)[[1]] %>% stringr::str_remove("\\|.*$")) & isFALSE(force)){
      stop("Accessions are  already in hexadecimal format, set force = TRUE to override this error")
      }
    acc <- names(x) %>%
      stringr::str_remove("\\|.*$") %>%
      purrr::map_chr(.hex)
    names(x) <- names(x) %>%
      stringr::str_remove("^.*(?=\\|)") %>%
      paste0(acc, .)
  }else  if (is(x, "character") && (stringr::str_detect(x, "\\|") & stringr::str_detect(x, ";"))) {
    message("Detected | and ; delimiters, assuming 'Accession|taxid;taxonomy' format")
    if(.ishex(x[[1]] %>% stringr::str_remove("\\|.*$")) & isFALSE(force)){
      stop("Accessions are  already in hexadecimal format, set force = TRUE to override this error")
    }
    acc <- x %>%
      stringr::str_remove("\\|.*$") %>%
      purrr::map_chr(.hex)
    x <- x %>%
      stringr::str_remove("^.*(?=\\|)") %>%
      paste0(acc, .)
  } else  if (is(x, "character") && !(stringr::str_detect(x, "\\|") && stringr::str_detect(x, ";"))) {
    message("Did not detect | and ; delimiters, assuming a vector of accessions")
    if(.ishex(x[[1]]) & isFALSE(force)){
      stop("Accessions are  already in hexadecimal format, set force = TRUE to override this error")
    }
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
#' @param force override checks if string is already hexadecimal
#' @return
#' @export
#' @import ape
#' @import Biostrings
#' @import stringr
#'
#' @examples
hex2acc <- function(x, force=FALSE){
  #Check input format
  if (is(x, "DNAbin") | is(x, "DNAStringSet")) {
    message("Input is ", class(x),", assuming header is in 'Accession|taxid;taxonomy' format")
    if(!.ishex(names(x)[[1]] %>% stringr::str_remove("\\|.*$")) & isFALSE(force)){
      stop("Accessions are already in alphanumeric, set force = TRUE to override this error")
    }
    acc <- names(x) %>%
      stringr::str_remove("\\|.*$") %>%
      purrr::map_chr(.unhex)
    names(x) <- names(x) %>%
      str_remove("^.*(?=\\|)") %>%
      paste0(acc, .)
  }else  if (is(x, "character") && (stringr::str_detect(x, "\\|") & stringr::str_detect(x, ";"))) {
    message("Detected | and ; delimiters, assuming 'Accession|taxid;taxonomy' format")
    if(!.ishex(x[[1]] %>% stringr::str_remove("\\|.*$")) & isFALSE(force)){
      stop("Accessions are already in alphanumeric, set force = TRUE to override this error")
    }
    acc <- x %>%
      stringr::str_remove("\\|.*$") %>%
      purrr::map_chr(.unhex)
    x <- x %>%
      stringr::str_remove("^.*(?=\\|)") %>%
      paste0(acc, .)
  } else  if (is(x, "character") && !(stringr::str_detect(x, "\\|") && stringr::str_detect(x, ";"))) {
    message("Did not detect | and ; delimiters, assuming a vector of accessions")
    if(!.ishex(x[[1]]) & isFALSE(force)){
      stop("Accessions are already in alphanumeric, set force = TRUE to override this error")
    }
    x <- x %>%
      purrr::map_chr(.unhex)
  }else (stop("x must be DNA bin, DNAStringSet or vector of accesssions"))
  return(x)
}


# hexadecimal translation -------------------------------------------------

#' convert alphanumerics to hexadecimal
#' internal taxreturn function
#' @param y a character string
#'
#' @return
#'
#' @examples
.hex <- function(y){
  paste(charToRaw(y),  collapse="")
}

#' convert hexadecimal strings to alphanumeric
#' internal taxreturn function
#' @param y a character string
#'
#' @return
#'
#' @examples
.unhex <- function(y){
  h <- sapply(seq(1, nchar(y), by=2), function(x) substr(y, x, x+1))
  rawToChar(as.raw(strtoi(h, 16L)))
}


#' Check whether a string is hexadecimal
#' internal taxreturn function
#' @param y a character string
#'
#' @return
#'
#' @examples
.ishex <- function(y) {
  h <- sapply(seq(1, nchar(y), by = 2), function(x) substr(y, x, x + 1))
  if(any(is.na(strtoi(h, 16L)))){
    return(FALSE)
  } else{
    return(TRUE)
  }
}
