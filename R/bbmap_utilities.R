#' Install BBtools
#'
#' @param url (Optional) Default will search for the latest version
#' URL to retrieve bbmap from.
#' @param dest.dir (Optional)  Default "bin"
#' Directory to install bbmap within.
#' @force Whether existing installs should be forcefully overwritten
#'
#'
#' @return
#' @export
#'
#' @import httr
#' @examples
bbmap_install <- function(url, dest.dir = "bin", force = FALSE) {
  if (missing(url)) {

    url <- ("https://sourceforge.net/projects/bbmap/files/latest/download")
  }

  if (!dir.exists(dest.dir)) {
    dir.create(dest.dir) # Create first directory
  }


  if (dir.exists(paste0(dest.dir, "/bbmap")) && force == FALSE) {
    stop("Stopped as bbmap already exists in directory, to overwrite set force to TRUE")
  } else  if (dir.exists(paste0(dest.dir, "/bbmap")) && force == TRUE) {
    unlink(paste0(dest.dir, "/bbmap"), recursive = TRUE) # Remove old version
  }

  destfile <- paste0(file.path(dest.dir, basename(url)),".tar.gz")
  if (file.exists(destfile)) {
    file.remove(destfile) # Remove old zip file
  }
  #Download file
  httr::GET(url, httr::write_disk(destfile, overwrite=TRUE))
  #Unzip
  utils::untar(destfile, exdir = dest.dir)
  #Remove download
  file.remove(destfile)
}


# Demultiplex by primers --------------------------------------------------

#' Demultiplex fusion primers using BBmap Seal
#'
#' @param install (Required) Install location for bbmap
#' @param fwd (Required) Vector of locations of forward reads
#' @param rev (Optional) Vector of locations of reverse reads
#' @param Fbarcodes (Required) Barcodes used in forward reads
#' @param Rbarcodes (Optional) Barcodes used in reverse reads
#' @param restrictleft (Optional) Defaults to the size of the largest primer.
#' Restricts the kmer search for primer sequences to just the left side of the molecule.
#' @param out.dir (Optional) Default "demux"
#'  The path to write the output reads.
#' @param kmer (Optional) default the size of the smallest primer will be used.
#' The kmer size to use for primer searching.
#' @param hdist (Optional) Default = 0. The hamming distance (number of substitution errors) allowed for mismatch to the query primer.
#' @param degenerate (Optional) Default TRUE.
#' Option to search for all possible primer combinations for degenerate primers
#' @param overwrite (Optional) Default TRUE
#' Option to overwrite existing output files.
#' @param interleaved (Optional) Default FALSE
#' Option to input interleaved reads
#' @param threads (Optional) Default autodetect
#' Number of CPU threads to use
#' WARNING: Thread detection can fail on cluster computing currently
#' @param mem (Optional) Default autodetect
#' GB of memory to use
#' WARNING: mem detection can fail on cluster computing currently
#'
#' @return
#' @export
#'
#' @import dplyr
#' @importFrom stringr str_split_fixed
#' @importFrom stringr str_replace
#'
#' @examples
#' #' \dontrun{
#' path <- "run_test/"
#' demuxpath <- file.path(path, "demux") # Filtered forward files go into the path/filtered/ subdirectory
#'
#' fastqFs <- sort(list.files(path, pattern="R1_001.*", full.names = TRUE))
#' fastqRs <- sort(list.files(path, pattern="R2_001.*", full.names = TRUE))
#'
#' bbdemux(install="bin/bbmap", fwd=fastqFs, rev=fastqRs,Fbarcodes = c("GAGGDACW","TGTGGDAC","AGAAGGDAC"),
#'               Rbarcodes = c("ACGTRATW","TCCGTRAT","CTGCGTRA"),
#'               degenerate=TRUE, out.dir=demuxpath, threads=1, mem=4,
#'               hdist=0, overwrite=TRUE)
#'               }
#'
bbdemux <- function(install = NULL, fwd, rev = NULL, Fbarcodes = NULL, Rbarcodes = NULL,
                          restrictleft = NULL, out.dir = "demux", kmer = NULL, hdist = 0, degenerate = TRUE,
                          overwrite = TRUE, threads = NULL, mem = NULL, interleaved = FALSE) {
  nsamples <- length(fwd)

  bbtools_seal <- function(install = NULL, fwd, rev = NULL, Fbarcodes = NULL, Rbarcodes = NULL,
                           restrictleft = NULL, out.dir = "demux", kmer = NULL, hdist = 0, degenerate = TRUE,
                           overwrite = TRUE, mem = NULL, threads = NULL) {
    in1 <- paste0("in=", fwd)
    if (!is.null(rev)) {
      in2 <- paste0("in2=", rev)
    } else {
      (in2 <- "")
    }

    if (!is.null(Fbarcodes) & is.null(Rbarcodes)) {
      writeLines(paste0(">Rep", seq(1:length(Fbarcodes)), "\n", Fbarcodes, "\n"), con = "Fprimers.fa")
      ref <- "ref=Fprimers.fa"
    } else if (!is.null(Fbarcodes) & !is.null(Rbarcodes)) {
      writeLines(paste0(">Rep", seq(1:length(Fbarcodes)), "\n", Fbarcodes, "\n"), con = "Fprimers.fa", sep = "")
      writeLines(paste0(">Rep", seq(1:length(Rbarcodes)), "\n", Rbarcodes, "\n"), con = "Rprimers.fa", sep = "")
      ref <- "ref=Fprimers.fa,Rprimers.fa"
    }

    pattern <- paste0("pattern=", out.dir, "/", basename(fwd) %>%
                        stringr::str_split_fixed("\\.", n = 2) %>%
                        tibble::as_tibble() %>%
                        dplyr::pull(V1) %>%
                        stringr::str_replace(pattern = "_R1_", replacement = "_R1R2_"), "_%.fastq.gz")

    if (is.numeric(kmer)) {
      kmer <- paste0("k=", kmer)
    } else {
      kmer <- paste0("k=", sort(nchar(c(Fbarcodes, Rbarcodes)), decreasing = FALSE)[1])
    }

    if (is.numeric(restrictleft)) {
      restrictleft <- paste0("restrictleft=", restrictleft)
    } else {
      restrictleft <- paste0("restrictleft=", sort(nchar(c(Fbarcodes, Rbarcodes)), decreasing = TRUE)[1])
    }


    if (is.numeric(hdist)) {
      hdist <- paste0("hdist=", hdist)
    }
    if (degenerate == TRUE) {
      degenerate <- "copyundefined"
    } else {
      (degenerate <- "")
    }

    if (overwrite == TRUE) {
      overwrite <- "overwrite=TRUE"
    } else {
      (overwrite <- "")
    }

    if (!is.null(threads)) {
      threads <- paste0("threads=", threads)
    } else {
      (threads <- "threads=auto")
    }

    if (!is.null(mem)) {
      mem <- paste0("-Xmx", mem, "g")
    } else {
      (mem <- "")
    }

    args <- paste(" -cp ", paste0(install, "/current jgi.Seal"), mem, in1, in2, ref,
                  restrictleft, pattern, kmer, hdist,
                  degenerate, overwrite, threads, "kpt=t",
                  collapse = " "
    )

    # Create temp files
    tmp <- tempdir()
    tmplogs <- paste0(tmp, "/bbdemux.log")
    tmpout <- paste0(tmp,"/stdout.log")
    tmperr <- paste0(tmp,"/stderr.log")


    result <- system2(command="java",
                      args = args,
                      stdout= tmpout,
                      stderr= tmperr,
                      wait=TRUE)
    now <- date()
    cat(paste0("Executed: ", now, "\n"), file=tmplogs, append=TRUE)
    cat(paste0("Sample:\t", fwd, "\n"), file=tmplogs, append=TRUE)
    file.append(tmplogs, tmperr)
    file.remove(c(tmpout, tmperr))
  }


  if (nsamples > 1) {
    for (i in 1:nsamples) {
      bbtools_seal(
        install = install, fwd = fwd[i], rev = rev[i], Fbarcodes = Fbarcodes, Rbarcodes = Rbarcodes,
        restrictleft = restrictleft, out.dir = out.dir, kmer = kmer, hdist = hdist, degenerate = degenerate,
        overwrite = overwrite, threads = threads, mem = mem
      )
    }
  } else if (nsamples == 1) {
    bbtools_seal(install, fwd, rev,
                 Fbarcodes = Fbarcodes, Rbarcodes = Rbarcodes, restrictleft = restrictleft,
                 out.dir = out.dir, kmer = kmer, hdist = hdist, degenerate = degenerate,
                 overwrite = overwrite, threads = threads, mem = mem
    )
  }
  #clean up
  file.remove("Fprimers.fa")
  file.remove("Rprimers.fa")

  #parse and return logs
  parsed <- parse_bbdemux(tmplogs)
  return(parsed)
}


# Trim primers ------------------------------------------------------------

# Note - may be worth adding an automated Maxlength to remove any untrimmed reads - maxlength = Readlength - shorterst primer + 1
# Note - Write out stats and output to a file within the output directory

#' Trim primers using BBDuk
#'
#' @param install (Required) Install location for bbmap
#' @param fwd (Required) Vector of locations of forward reads
#' @param rev (Optional) Vector of locations of reverse reads
#' @param primers (Required) Forward and reverse primers to trim
#' @param restrictleft (Optional) Defaults to the size of the largest primer.
#' Restricts the kmer search for primer sequences to just the left side of the molecule.
#' @param out.dir (Optional) Default "trimmed"
#'  The path to write the output reads.
#' @param trim.end (Optional) Default is "left"
#' End of the molecule to trim primers from. "left" will trim primers from the 3' end of both forward and reverse reads.
#' Change to "right" only if the amplicon was too short and the sequencer has read into the other end of the molecule.
#' @param ordered (Optional) Default TRUE
#'  Set to TRUE to output reads in same order as input.
#' @param kmer (Optional) default the size of the smallest primer will be used.
#' The kmer size to use for primer searching.
#' @param mink (optional) Default FALSE
#' Look for shorter kmers at read tips down to this length
#' @param hdist (Optional) Default 0. The hamming distance (number of substitution errors) allowed for mismatch to the query primer.
#' @param tpe (Otional) Default TRUE
#' Trim pairs evenly. When kmer right-trimming, trim both reads to the minimum length of either.
#' @param degenerate (Optional) Default TRUE.
#' Option to search for all possible primer combinations for degenerate primers
#' @param overwrite (Optional) Default TRUE
#' Option to overwrite existing output files.
#' @param quality (Optional) Default FALSE
#' Output quality statistics from trimming including:
#' Base composition histogram by position, Quality histogram by position,
#' Count of bases with each quality value, Histogram of average read quality,
#' Read length histogram and Read GC content histogram.
#' @param maxlength (Optional) Default FALSE
#' Remove all reads above a maximum length. Useful for removing reads where no primers were found.
#'
#' @return
#' @export
#'
#' @import dplyr
#' @import tidyr
#' @import readr
#' @import purrr
#' @importFrom stringr str_replace
#'
#'
#' @examples
#' \dontrun{
#' path <- "run_test/"
#'
#' fastqFs <- sort(list.files(path, pattern="R1_001.*", full.names = TRUE))
#' fastqRs <- sort(list.files(path, pattern="R2_001.*", full.names = TRUE))
#'
#'bbtrim(install="bin/bbmap", fwd=fastqFs, rev=fastqRs,
#' primers=c("GGDACWGGWTGAACWGTWTAYCCHCC","GTRATWGCHCCDGCTARWACWGG"),
#'  degenerate=TRUE, out.dir="trimmed", ktrim="left", ordered=TRUE,
#'   mink=FALSE, hdist=2, maxlength=140, overwrite=TRUE)
#'}
bbtrim <- function(install = NULL, fwd, rev = NULL, primers,
                         restrictleft = NULL, out.dir = "bbduk", trim.end = "left", ordered = TRUE,
                         kmer = NULL, mink = FALSE, tpe = TRUE, hdist = 0, degenerate = TRUE,
                         overwrite = TRUE, quality = FALSE, maxlength = NULL) {
  nsamples <- length(fwd)

      bbduk <- function(install = NULL, fwd, rev = NULL, primers,
                        restrictleft = NULL, out.dir = "bbduk", trim.end = "left", ordered = TRUE,
                        kmer = NULL, mink = FALSE, tpe = TRUE, hdist = 0, degenerate = TRUE,
                        overwrite = TRUE, quality = FALSE, maxlength = NULL) {
        install <- paste0(install, "/current jgi.BBDuk")

        in1 <- paste0("in=", fwd)
        if (!is.null(rev)) {
          in2 <- paste0("in2=", rev)
        } else {
          (in2 <- "")
        }

        if (!is.null(primers)) {
          literal <- paste0("literal=", paste0(primers, collapse = ","))
        } else {
          (stop("Primer sequences are required for trimming"))
        }

        if (is.null(rev)) {
          out <- paste0(
            "out=", dirname(fwd), "/", out.dir, "/", basename(fwd)
          ) %>%
            stringr::str_replace(pattern = ".fastq", replacement = ".trimmed.fastq")
          out1 <- ""
          out2 <- ""
        } else if (!is.null(rev)) {
          out <- ""
          out1 <- paste0("out1=", dirname(fwd), "/", out.dir, "/", basename(fwd)) %>%
            stringr::str_replace(pattern = ".fastq", replacement = ".trimmed.fastq")
          out2 <- paste0("out2=", dirname(rev), "/", out.dir, "/", basename(rev)) %>%
            stringr::str_replace(pattern = ".fastq", replacement = ".trimmed.fastq")
        }


        if (trim.end == "left") {
          trim.end <- paste0("ktrim=l")
        } else if (trim.end == "right") {
          trim.end <- paste0("ktrim=r")
        }

        if (is.numeric(kmer)) {
          kmer <- paste0("k=", kmer)
        } else {
          (kmer <- paste0("k=", sort(nchar(primers), decreasing = FALSE)[1]))
        }

        if (is.numeric(maxlength)) {
          maxlength <- paste0("maxlength=", maxlength)
        } else {
          maxlength <- ""
        }

        if (is.numeric(mink)) {
          mink <- paste0("mink=", mink) # Note - mink makes it noticibly slower
        } else if (mink == TRUE) {
          mink <- paste0("mink=", (sort(nchar(primers), decreasing = FALSE)[1] / 2))
        } else if (mink == FALSE) {
          mink <- ""
        }


        if (is.numeric(restrictleft)) {
          restrictleft <- paste0("restrictleft=", restrictleft)
        } else {
          restrictleft <- paste0("restrictleft=", sort(nchar(primers), decreasing = TRUE)[1])
        }

        if (ordered == TRUE) {
          ordered <- "ordered=T"
        } else {
          (ordered <- "")
        }

        if (is.numeric(hdist)) {
          hdist <- paste0("hdist=", hdist)
        }

        if (degenerate == TRUE) {
          degenerate <- "copyundefined"
        } else {
          (degenerate <- "")
        }

        if (overwrite == TRUE) {
          overwrite <- "overwrite=TRUE"
        } else {
          (overwrite <- "")
        }

        if (tpe == TRUE) {
          tpe <- "tpe"
        } else {
          (tpe <- "")
        }

        args <- paste(" -cp ", install, in1, in2, literal, restrictleft, out, out1,
          out2, kmer, mink, hdist, trim.end, tpe, degenerate, quality,
          maxlength, overwrite, "-da",
          collapse = " "
        )

        # Create temp files
        tmp <- tempdir()
        tmplogs <- paste0(tmp, "/bbtrim.log")
        tmpout <- paste0(tmp,"/stdout.log")
        tmperr <- paste0(tmp,"/stderr.log")

        # Set up quality tracking
        if (quality == TRUE) {
          qualnames <- fwd %>% str_replace(pattern=".fastq.gz", replacement="") %>%
            basename
          qualnames <- paste0(tmp, qualnames)
          quality <-
            paste0("bhist=", qualnames, "_bhist.txt ",
                   "qhist=", qualnames, "_qhist.txt ",
                   "gchist=", qualnames, "_gchist.txt ",
                   "aqhist=", qualnames, "_aqhist.txt ",
                   "lhist=", qualnames, "_lhist.txt ",
                   "gcbins=auto ")
        } else {
          (quality <- "")
        }

        # Run bbduk
        result <- system2(command="java",
                          args = args,
                          stdout = tmpout,
                          stderr = tmperr,
                          wait=TRUE)
        now <- date()
        cat(paste0("Executed: ", now, "\n"), file = tmplogs, append=TRUE)
        cat(paste0("Sample:\t", fwd, "\n"), file = tmplogs, append=TRUE)
        file.append(tmplogs, tmperr)
        file.remove(c(tmpout, tmperr))
      }

  if (nsamples > 1) {
    for (i in 1:nsamples) {
      bbduk(install = install, fwd = fwd[i], rev = rev[i],
            primers = primers, restrictleft = restrictleft,
            out.dir = out.dir, trim.end = trim.end, ordered = ordered,
            kmer = kmer, mink = mink, tpe = tpe, hdist = hdist,
            degenerate = degenerate, quality = quality,
            overwrite = overwrite, maxlength = maxlength)
    }
  } else if (nsamples == 1) {
    bbduk(install = install, fwd = fwd, rev = rev,
          primers = primers, restrictleft = restrictleft,
          out.dir = out.dir, trim.end = trim.end, ordered = ordered,
          kmer = kmer, mink = mink, tpe = tpe, hdist = hdist,
          degenerate = degenerate, quality = quality,
          overwrite = overwrite, maxlength = maxlength)
  }

  #Parse logs
  parsed <- parse_bbtrim(tmplogs)

  if (quality == TRUE) {
    #Base composition histogram by position.
    bhist <- parse_bhist(tmp)
    #Quality histogram by position.
    qhist <- parse_qhist(tmp)
    #Histogram of average read quality. - how does this work with binned qscores?
    aqhist <- parse_aqhist(tmp)
    #Read GC content histogram. - is it worth just reading in the top 4 lines?
    gchist <- parse_aqhist(tmp)
    #Read length histogram.
    lhist <- parse_lhist(tmp)

    out <- list(parsed,
                bhist,
                qhist,
                aqhist,
                gchist,
                lhist)

  } else (out <- parsed)

    return(out)
}

# Split interleaved reads -------------------------------------------------

#' Split interleaved reads
#'
#' @param install (Required) Install location for bbmap
#'
#' @param files (Required) Vector of locations of interleaved read files to split
#'
#' @param overwrite (Optional) Default TRUE
#' Option to overwrite existing output files.
#'
#' @return
#' @export
#'
bbsplit <- function(install = NULL, files, overwrite = FALSE) {
  nsamples <- length(files)

  bbtools_reformatreads <- function(install = NULL, file, overwrite = FALSE) {
    # Split interleaved reads
    out1 <- paste0("out1=", file %>% str_replace(pattern = "_R1R2_", replacement = "_R1_"))
    out2 <- paste0("out2=", file %>% str_replace(pattern = "_R1R2_", replacement = "_R2_"))

    if (overwrite == TRUE) {
      overwrite <- "overwrite=TRUE"
    } else {
      (overwrite <- "")
    }

    # Create temp files
    tmp <- tempdir()
    tmplogs <- paste0(tmp, "/bbreformat.log")
    tmpout <- paste0(tmp,"/stdout.log")
    tmperr <- paste0(tmp,"/stderr.log")

    reformat_args <- paste(" -cp ", paste0(install, "/current jgi.ReformatReads "), file, out1, out2, overwrite, collapse = " ")

    # Run Reformatreads
    result <- system2(command="java",
                      args = reformat_args,
                      stdout=tmpout,
                      stderr=tmperr,
                      wait=TRUE)
    now <- date()
    cat(paste0("Executed: ", now, "\n"), file="logs/bbreformat.log", append=TRUE)
    file.append(tmplogs, tmperr)
    file.remove(c(tmpout, tmperr))

  }

  if (nsamples > 1) {
    for (i in 1:nsamples) {
      bbtools_reformatreads(install = install, file = files[i], overwrite = overwrite)
      file.remove(files[i])
    }
  } else if (nsamples == 1) {
    bbtools_reformatreads(install = install, file = files, overwrite = overwrite)
    file.remove(files)
  }
}



# Parse bbtrim logs ----------------------------------------------------------

#' Parse appended bbtrim logs into tidy format
#'
#' @param x (Required) an appended bbtrim log file generated by bbtrim
#'
#' @return Returns a tidy data frame
#'
#' @export
#'
#' @examples
parse_bbtrim <- function(x) {

  lines <- readLines(x)

  sample <- readr::read_tsv(lines[which(str_sub(lines, 1, 7) == 'Sample:' )], col_names = FALSE) %>%
    tidyr::separate(X2, into="sample", sep="\\.", extra="drop") %>%
    dplyr::select(sample)

  input <- read_tsv(lines[which(str_sub(lines, 1, 6) == 'Input:' )], col_names = FALSE) %>%
    tidyr::separate(X2, into="input_reads", sep=" ", extra="drop") %>%
    tidyr::separate(X4, into="input_bases", sep=" ", extra="drop") %>%
    dplyr::select(input_reads, input_bases) %>%
    dplyr::mutate_if(is.character, as.numeric)

  ktrimmed <- read_tsv(lines[which(str_sub(lines, 1, 9) == 'KTrimmed:' )], col_names = FALSE) %>%
    tidyr::separate(X2, into="ktrimmed_reads", sep=" ", extra="drop") %>%
    tidyr::separate(X3, into="ktrimmed_bases", sep=" ", extra="drop") %>%
    dplyr::select(ktrimmed_reads, ktrimmed_bases) %>%
    dplyr::mutate_if(is.character, as.numeric)

  result <- read_tsv(lines[which(str_sub(lines, 1, 7) == 'Result:' )], col_names = FALSE) %>%
    tidyr::separate(X2, into="output_reads", sep=" ", extra="drop") %>%
    tidyr::separate(X3, into="output_bases", sep=" ", extra="drop") %>%
    dplyr::select(reads_out, bases_out) %>%
    dplyr::mutate_if(is.character, as.numeric)

  out <- cbind(sample, input, ktrimmed, result)
  return(out)

}

# Parse bbdemux -----------------------------------------------------------

#' Parse appended bbdemux logs into tidy format
#'
#' @param x (Required) an appended bbdemux log file generated by bbdemux
#'
#' @return Returns a tidy data frame
#'
#' @examples
parse_bbdemux <- function(x) {

  lines <- readLines(x)

  sample <- readr::read_tsv(lines[which(str_sub(lines, 1, 7) == 'Sample:' )], col_names = FALSE) %>%
    tidyr::separate(X2, into="sample", sep="\\.", extra="drop") %>%
    dplyr::select(sample)

  input <- readr::read_tsv(lines[which(str_sub(lines, 1, 6) == 'Input:' )], col_names = FALSE) %>%
    tidyr::separate(X2, into="input_reads", sep=" ", extra="drop") %>%
    tidyr::separate(X4, into="input_bases", sep=" ", extra="drop") %>%
    dplyr::select(input_reads, input_bases) %>%
    dplyr::mutate_if(is.character, as.numeric)

  matched <- readr::read_tsv(lines[which(str_sub(lines, 1, 14) == 'Matched reads:' )], col_names = FALSE) %>%
    tidyr::separate(X2, into="demulti_reads", sep=" ", extra="drop") %>%
    tidyr::separate(X3, into="demulti_bases", sep=" ", extra="drop") %>%
    dplyr::select(demulti_reads, demulti_bases) %>%
    dplyr::mutate_if(is.character, as.numeric)

 out <-  cbind(sample, input, matched )

 return(out)
}


# parse_bhist -------------------------------------------------------------

#' parse base frequency histograms from bbtrim
#'
#' @param dir (Required) a directory containing base frequency (bhist) files from bbtrim
#'
#' @return a tidy data frame
#'
#' @examples
parse_bhist <- function(dir) {
  out <- list.files(path=dir, pattern="_bhist.txt", full.names = TRUE) %>%
    purrr::set_names()  %>%
    purrr::map_dfr(readr::read_table2, .id = "Source", col_types = cols(
      `#Pos` = col_double(),
      A = col_double(),
      C = col_double(),
      G = col_double(),
      T = col_double(),
      N = col_double()
    )) %>%
    dplyr::rename(Pos = `#Pos`)

  return(out)
}


# parse_qhist -------------------------------------------------------------
#' Parse quality score histograms from bbtrim
#'
#' @param dir (Required) a directory containing quality score (qhist) files from bbtrim
#'
#' @return a tidy data frame
#'
#' @examples
parse_qhist <- function(dir) {
  out  <- list.files(path=dir, pattern="_qhist.txt", full.names = TRUE) %>%
    purrr::set_names()  %>%
    purrr::map_dfr(readr::read_table2, .id = "Source", col_types = cols(
      `#BaseNum` = col_double(),
      Read1_linear = col_double(),
      Read1_log = col_double(),
      Read2_linear = col_double(),
      Read2_log = col_double()
    )) %>%
    dplyr::rename(Pos = `#BaseNum`)

  return(out)
}


# parse_aqhist ------------------------------------------------------------
#' Parse average quality histograms from bbtrim
#'
#' @param dir (Required) a directory containing average quality (aqhist) files from bbtrim
#'
#' @return a tidy data frame
#'
#' @examples
parse_aqhist <- function(dir) {
  out  <- list.files(path=dir, pattern="_aqhist.txt", full.names = TRUE) %>%
    purrr::set_names()  %>%
    purrr::map_dfr(readr::read_table2, .id = "Source", col_types = cols(
      `#Quality` = col_double(),
      count1 = col_double(),
      fraction1 = col_double(),
      count2 = col_double(),
      fraction2 = col_double()
    ))%>%
    dplyr::rename(avg_quality = `#Quality`)

  return(out)
}


# parse_gchist ------------------------------------------------------------

#' Parse GC content histograms from bbtrim
#'
#' @param dir (Required) a directory containing average quality (aqhist) files from bbtrim
#'
#' @return a tidy data frame
#'
#' @examples
parse_gchist <- function(dir) {
  out <- list.files(path=dir, pattern="_gchist.txt", full.names = TRUE) %>%
    purrr::set_names()  %>%
    purrr::map_dfr(readr::read_table2, .id = "Source", n_max=4, col_names = FALSE, col_types = cols(
      X1 = col_character(),
      X2 = col_double()
    ))  %>%
    dplyr::mutate(X1 = stringr::str_replace(`X1`, pattern="#", replacement="")) %>%
    dplyr::group_by(Source) %>%
    tidyr::pivot_wider(names_from = X1, values_from = X2)

  return(out)
}


# parse_lhist -------------------------------------------------------------
#' Parse sequence length histograms from bbtrim
#'
#' @param dir (Required) a directory containing length (lhist) files from bbtrim
#'
#' @return a tidy data frame
#'
#' @examples
parse_lhist <- function(dir) {
  out  <- list.files(path=dir, pattern="_lhist.txt", full.names = TRUE) %>%
    purrr::set_names()  %>%
    purrr::map_dfr(readr::read_table2, .id = "Source", col_types= cols(
      `#Length` = col_double(),
      Count = col_double()
    )) %>%
    dplyr::rename(Length = `#Length`)
}
