# DNAbin2DNAStringSet -----------------------------------------------------

#' Convert DNABin to DNAStringSet
#'
#' @param x a DNABin object
#' @param remove_gaps Whether gaps should be removed
#'
#' @return
#' @export
#' @import purrr
#' @importFrom Biostrings DNA_ALPHABET
#' @importFrom Biostrings DNAStringSet
#' @importFrom methods is
#'
DNAbin2DNAstringset <- function (x, remove_gaps = FALSE) {
  if(!methods::is(x, "DNAbin")){
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
#' @importFrom insect char2dna
#'
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
#' @import stringr
#' @importFrom methods is
#'
acc2hex <- function(x, force=FALSE){
  #Check input format
  if (methods::is(x, "DNAbin") | methods::is(x, "DNAStringSet")) {
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
  }else  if (methods::is(x, "character") && (stringr::str_detect(x, "\\|") & stringr::str_detect(x, ";"))) {
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
  } else  if (methods::is(x, "character") && !(stringr::str_detect(x, "\\|") && stringr::str_detect(x, ";"))) {
    message("Did not detect | and ; delimiters, assuming a vector of accessions")
    if(.ishex(x[[1]]) & isFALSE(force)){
      stop("Accessions are  already in hexadecimal format, set force = TRUE to override this error")
    }
    x <- x %>%
      purrr::map_chr(.hex)
  }else (stop("x must be DNA bin, DNAStringSet or vector of accesssions"))
  return(x)
}

#' Concatenate DNAbin objects while preserving attributes.
#' This function joins two or more \code{DNAbin} objects, retaining any
#'   attributes whose lengths match those of the input objects (e.g. "species",
#'   "lineage" and/or "taxID" attributes).
#' @param ... \code{DNAbin} objects, or list of DNAbins to be concatenated.
#'
#' @return
#' @export
#'
#' @examples
concat_DNAbin <- function(...){
  dots <- list(...)
  if(is.list(dots[[1]]) & length(dots) == 1){
    dots <- dots[[1]]
  }
  to_remove <- sapply(dots, is.null)
  dots <- dots[!to_remove]
  nlsts <- length(dots)
  DNA <- any(sapply(dots, class) == "DNAbin")
  if(nlsts == 0) return(NULL)
  if(nlsts == 1) return(dots[[1]])
  islist <- sapply(dots, is.list)
  for(i in which(!islist)){
    print(i)
    tmpattr <- attributes(dots[[i]])
    attributes(dots[[i]]) <- NULL
    dots[[i]] <- list(dots[[i]])

    print(tmpattr)
    attributes(dots[[i]]) <- tmpattr
  }
  dots <- lapply(dots, unclass)
  findattr <- function(x){
    names(attributes(x))[sapply(attributes(x), length) == length(x)]
  }
  attrlist <- lapply(dots, findattr)
  ual <- unique(unlist(attrlist, use.names = FALSE))
  validattrs <- ual[sapply(ual, function(e) all(sapply(attrlist, function(g) e %in% g)))]
  validattrs <- validattrs[!validattrs %in% c("names", "class")]
  res <- unlist(dots, recursive = FALSE, use.names = TRUE)
  for(i in validattrs){
    attr(res, i) <- unlist(lapply(dots, attr, i), use.names = FALSE)
  }
  class(res) <- "DNAbin"
  return(res)
}


# Hexadecimal coding to accession -----------------------------------------
#' Hexadecimal coding to accession numberD
#' @description Hexadecimal coding is useful for ensuring that accession numbers dont contain characters that are considered illegal for newick strings
#' @param x A DNAbin or DNAStringSet with names formatted Accession|taxid;taxonomy.
#' Or a character vector with names formatted Accession|taxid;taxonomy
#' Or a character vector of accessions
#' @param force override checks if string is already hexadecimal
#' @return
#' @export
#' @import stringr
#' @importFrom methods is
#'
hex2acc <- function(x, force=FALSE){
  #Check input format
  if (methods::is(x, "DNAbin") | methods::is(x, "DNAStringSet")) {
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
  }else  if (methods::is(x, "character") && (stringr::str_detect(x, "\\|") & stringr::str_detect(x, ";"))) {
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
  } else  if (methods::is(x, "character") && !(stringr::str_detect(x, "\\|") && stringr::str_detect(x, ";"))) {
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
#'
.hex <- function(y){
  paste(charToRaw(y),  collapse="")
}

#' convert hexadecimal strings to alphanumeric
#' internal taxreturn function
#' @param y a character string
#'
#' @return
#'
#'
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
#'
.ishex <- function(y) {
  h <- sapply(seq(1, nchar(y), by = 2), function(x) substr(y, x, x + 1))
  if(any(is.na(strtoi(h, 16L)))){
    return(FALSE)
  } else{
    return(TRUE)
  }
}



# Multithread -------------------------------------------------------------
#' Setup multithreading
#'
#' @param multithread Number of cores
#' @param quiet Whether progress should be printed to console
#'
#' @return
#' @import future
#'
#' @examples
setup_multithread <- function(multithread, quiet=FALSE){
  ncores <- future::availableCores()
  if(isTRUE(multithread)){
    cores <- ncores-1
    if(!quiet){message("Multithreading with ", cores, " cores")}
    future::plan(future::multiprocess, workers=cores)
  } else if (is.numeric(multithread) & multithread > 1){
    cores <- multithread
    if(cores > ncores){
      cores <- ncores
      warning("The value provided to multithread is higher than the number of cores, using ", cores, " cores instead")
    }
    if(!quiet){message("Multithreading with ", cores, " cores")}
    future::plan(future::multiprocess, workers=cores)
  } else if(isFALSE(multithread) | multithread==1){
    future::plan(future::sequential)
  } else (
    stop("Multithread must be a logical or numeric vector of the numbers of cores to use")
  )
}

#' Setup parallel
#'
#' @param multithread Number of cores
#' @param quiet Whether progress should be printed to console
#'
#' @return
#' @importFrom future availableCores
#'
#' @examples
setup_para <- function(multithread, quiet=FALSE){
  ncores <- future::availableCores() -1
  if(isTRUE(multithread)){
    cores <- ncores
    if(!quiet){message("Multithreading with ", cores, " cores")}
  } else if (is.numeric(multithread) & multithread > 1){
    cores <- multithread
    if(cores > ncores){
      cores <- ncores
      warning("The value provided to multithread is higher than the number of cores, using ", cores, " cores instead")
    }
    if(!quiet){message("Multithreading with ", cores, " cores")}
  } else if(isFALSE(multithread) | multithread==1){
    cores <- 1
  } else (
    stop("Multithread must be a logical or numeric vector of the numbers of cores to use")
  )
  return(cores)
}
