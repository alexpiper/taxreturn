# DNAbin Utilities -----------------------------------------------------

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

#' char2DNAbin
#'
#' @param x a character vector
#'
#' @return
#' @export
#'
#' @examples
char2DNAbin <- function (x) {
  dbytes <- as.raw(c(136, 24, 72, 40, 96, 144, 192, 48, 80,
                     160, 112, 224, 176, 208, 240, 240, 4, 2))
  indices <- c(65, 84, 71, 67, 83, 87, 82, 89, 75, 77, 66,
               86, 72, 68, 78, 73, 45, 63)
  vec <- raw(89)
  vec[indices] <- dbytes
  s2d1 <- function(s) vec[as.integer(charToRaw(s))]
  out <- lapply(x, s2d1)
  attr(out, "rerep.names") <- attr(x, "rerep.names")
  attr(out, "rerep.pointers") <- attr(x, "rerep.pointers")
  class(out) <- "DNAbin"
  return(out)
}

#' DNAbin to char
#'
#' @param x A DNAbin object
#'
#' @return
#' @export
#'
#' @examples
DNAbin2char <- function(x){
  cbytes <- as.raw(c(65, 84, 71, 67, 83, 87, 82, 89, 75, 77,
                     66, 86, 72, 68, 78, 45, 63))
  indices <- c(136, 24, 72, 40, 96, 144, 192, 48, 80, 160,
               112, 224, 176, 208, 240, 4, 2)
  vec <- raw(240)
  vec[indices] <- cbytes
  if (is.list(x)) {
    out <- sapply(x, function(s) rawToChar(vec[as.integer(s)]))
  } else {
    out <- rawToChar(vec[as.integer(x)])
  }
  attr(out, "rerep.names") <- attr(x, "rerep.names")
  attr(out, "rerep.pointers") <- attr(x, "rerep.pointers")
  return(out)
}

# add a queit arg
# add a detect "gz" set compress to true

#' write fasta
#'
#' @param x a list of sequences in DNAbin or AAbin format, or a vector of sequences as concatenated upper-case character strings.
#' @param file character string giving a valid file path to output the text to. If file = "" (default setting) the text file is written to the console.
#' @param compress logical indicating whether the output file should be gzipped.
#' @param quiet Whether progress should be printed to consoe
#'
#' @return
#' @export
#'
#' @examples
write_fasta <- function(x, file = "", compress = FALSE, quiet=FALSE) {
  if(stringr::str_detect(file, "\\.gz$") & !compress){
    compress <- TRUE
    if(!quiet) message(".gz detected in filename, compressing output file")
  }
  if (!is.null(dim(x))) {
    x <- as.list(as.data.frame(t(unclass(x))))
  }
  if (inherits(x, "DNAbin")) {
    tmp <- DNAbin2char(x)
  } else if (is.list(x)) {
    if (length(x[[1]] == 1)) {
      tmp <- unlist(x, use.names = TRUE)
    }
    else {
      tmp <- sapply(x, paste0, collapse = "")
    }
  } else {
    tmp <- x
  }
  reslen <- 2 * length(tmp)
  res <- character(reslen)
  res[seq(1, reslen, by = 2)] <- paste0(">", names(tmp))
  res[seq(2, reslen, by = 2)] <- tmp
  if(!file == ""){
    f <- if(compress){
        gzfile(file, "w")
      } else {
        file(file, "w")
    }
    writeLines(res, f)
    close(f)
  } else {
    writeLines(res)
  }
  if(!quiet){message("Wrote ", length(tmp), " seqeuences to ", file)}
  invisible(NULL)
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
