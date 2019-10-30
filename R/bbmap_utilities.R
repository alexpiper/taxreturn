#' BBMap Install
#'
#' @param url (Optional) Default = "https://downloads.sourceforge.net/project/bbmap/BBMap_38.57.tar.gz"
#' URL to retrieve bbmap from.
#' @param destdir (Optional)  Default "bin"
#' Directory to install bbmap within.
#'
#' @return
#' @export
#'
#' @examples
bbmap_install <- function(url, destdir = "bin") {
  if (missing(url)) {
    url <- ("https://downloads.sourceforge.net/project/bbmap/BBMap_38.57.tar.gz")
  }

  if (!dir.exists(destdir)) {
    dir.create(destdir)
  } # Create first directory

  if (dir.exists(paste0(destdir, "/bbmap"))) {
    unlink(paste0(destdir, "/bbmap"), recursive = TRUE)
  } # Remove old version

  destfile <- file.path(destdir, basename(url))
  if (exists(destfile)) {
    file.remove(destfile)
  } # Remove old zip file

  utils::download.file(url, destfile = destfile)
  utils::untar(destfile, exdir = destdir) ## check contents
  file.remove(destfile)
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
#' @param outpath (Optional) Default "trimmed"
#'  The path to write the output reads.
#' @param trim.dir (Optional) Default is "left"
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
#' Output quality statistics from trimming
#' @param maxlength (Optional) Default FALSE
#' Remove all reads above a maximum length. Useful for removing reads where no primers were found.
#'
#' @return
#' @export
#'
#'@importFrom stringr str_replace
#'
#' @examples
#' path <- "run_test/"
#' demuxpath <- file.path(path, "demux") # Filtered forward files go into the path/filtered/ subdirectory
#'
#' fastqFs <- sort(list.files(path, pattern="R1_001.*", full.names = TRUE))
#' fastqRs <- sort(list.files(path, pattern="R2_001.*", full.names = TRUE))
#'
#'bbtools_trim(install="bin/bbmap", fwd=fastqFs, rev=fastqRs,
#' primers=c("GGDACWGGWTGAACWGTWTAYCCHCC","GTRATWGCHCCDGCTARWACWGG"),
#'  degenerate=TRUE, outpath="trimmed", ktrim="left", ordered=TRUE,
#'   mink=FALSE, hdist=2, maxlength=140, overwrite=TRUE)
#'
bbtools_trim <- function(install = NULL, fwd, rev = NULL, primers,
                         restrictleft = NULL, outpath = "bbduk", trim.dir = "left", ordered = TRUE,
                         kmer = NULL, mink = FALSE, tpe = TRUE, hdist = 0, degenerate = TRUE,
                         overwrite = TRUE, quality = FALSE, maxlength = NULL) {
  nsamples <- length(fwd)

  bbduk <- function(install = NULL, fwd, rev = NULL, primers,
                    restrictleft = NULL, outpath = "bbduk", trim.dir = "left", ordered = TRUE,
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
        "out=", dirname(fwd), "/", outpath, "/", basename(fwd)
      ) %>%
        stringr::str_replace(pattern = ".fastq", replacement = ".trimmed.fastq")
      out1 <- ""
      out2 <- ""
    } else if (!is.null(rev)) {
      out <- ""
      out1 <- paste0("out1=", dirname(fwd), "/", outpath, "/", basename(fwd)) %>%
        stringr::str_replace(pattern = ".fastq", replacement = ".trimmed.fastq")
      out2 <- paste0("out2=", dirname(rev), "/", outpath, "/", basename(rev)) %>%
        stringr::str_replace(pattern = ".fastq", replacement = ".trimmed.fastq")
    }


    if (trim.dir == "left") {
      trim.dir <- paste0("ktrim=l")
    } else if (trim.dir == "right") {
      trim.dir <- paste0("ktrim=r")
    }

    if (is.numeric(kmer)) {
      kmer <- paste0("k=", kmer)
    } else {
      kmer <- paste0("k=", nchar(sort(primers, decreasing = TRUE)[1]))
    }

    if (is.numeric(maxlength)) {
      maxlength <- paste0("maxlength=", maxlength)
    } else {
      maxlength <- ""
    }

    if (is.numeric(mink)) {
      mink <- paste0("mink=", mink) # Note - mink makes it noticibly slower
    } else if (mink == TRUE) {
      mink <- paste0("mink=", (ceiling(nchar(sort(primers, decreasing = TRUE)[1]) / 2)))
    } else if (mink == FALSE) {
      mink <- ""
    }


    if (is.numeric(restrictleft)) {
      restrictleft <- paste0("restrictleft=", restrictleft)
    } else {
      (restrictleft <- paste0("restrictleft=", nchar(sort(primers, decreasing = FALSE)[1]) + 1))
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

    if (quality == TRUE) {
      quality <- "bhist=bhist.txt qhist=qhist.txt gchist=gchist.txt aqhist=aqhist.txt lhist=lhist.txt gcbins=auto"
    } else {
      (quality <- "")
    }

    args <- paste(" -cp ", install, in1, in2, literal, restrictleft, out, out1,
      out2, kmer, mink, hdist, trim.dir, tpe, degenerate, quality,
      maxlength, overwrite, "-da >> stdout.log 2>> stderr.log",
      collapse = " "
    )

    print(paste0(args, " /n"))
    # Run bbduk
    system2("java", args = args, wait = TRUE)
  }
  if (nsamples > 1) {
    for (i in 1:nsamples) {
      bbduk(install = install, fwd = fwd[i], rev = rev[i],
            primers = primers, restrictleft = restrictleft,
            outpath = outpath, trim.dir = trim.dir, ordered = ordered,
            kmer = kmer, mink = mink, tpe = tpe, hdist = hdist,
            degenerate = degenerate, quality = quality,
            overwrite = overwrite, maxlength = maxlength)
    }
  } else if (nsamples == 1) {
    bbduk(install = install, fwd = fwd, rev = rev,
          primers = primers, restrictleft = restrictleft,
          outpath = outpath, trim.dir = trim.dir, ordered = ordered,
          kmer = kmer, mink = mink, tpe = tpe, hdist = hdist,
          degenerate = degenerate, quality = quality,
          overwrite = overwrite, maxlength = maxlength)
  }
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
#' @param outpath (Optional) Default "demux"
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
#' @importFrom stringr str_split_fixed
#' @importFrom  stringr str_replace
#' @importFrom tibble as_tibble
#' @importFrom dplyr pull
#'
#' @examples
#' path <- "run_test/"
#' demuxpath <- file.path(path, "demux") # Filtered forward files go into the path/filtered/ subdirectory
#'
#' fastqFs <- sort(list.files(path, pattern="R1_001.*", full.names = TRUE))
#' fastqRs <- sort(list.files(path, pattern="R2_001.*", full.names = TRUE))
#'
#' bbtools_demux(install="bin/bbmap", fwd=fastqFs, rev=fastqRs,Fbarcodes = c("GAGGDACW","TGTGGDAC","AGAAGGDAC"),
#'               Rbarcodes = c("ACGTRATW","TCCGTRAT","CTGCGTRA"),
#'               degenerate=TRUE, outpath=demuxpath, threads=1, mem=4,
#'               hdist=0, overwrite=TRUE)
#'
bbtools_demux <- function(install = NULL, fwd, rev = NULL, Fbarcodes = NULL, Rbarcodes = NULL,
                          restrictleft = NULL, outpath = "demux", kmer = NULL, hdist = 0, degenerate = TRUE,
                          overwrite = TRUE, threads = NULL, mem = NULL, interleaved = FALSE) {
  nsamples <- length(fwd)

  bbtools_seal <- function(install = NULL, fwd, rev = NULL, Fbarcodes = NULL, Rbarcodes = NULL,
                           restrictleft = NULL, outpath = "demux", kmer = NULL, hdist = 0, degenerate = TRUE,
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

    pattern <- paste0("pattern=", outpath, "/", basename(fwd) %>%
      stringr::str_split_fixed("\\.", n = 2) %>%
      tibble::as_tibble() %>%
      dplyr::pull(V1) %>%
      stringr::str_replace(pattern = "_R1_", replacement = "_R1R2_"), "_%.fastq.gz")

    if (is.numeric(kmer)) {
      kmer <- paste0("k=", kmer)
    } else {
      (kmer <- paste0("k=", nchar(sort(c(Fbarcodes, Rbarcodes), decreasing = TRUE)[1])))
    }

    if (is.numeric(restrictleft)) {
      restrictleft <- paste0("restrictleft=", restrictleft)
    } else {
      (restrictleft <- paste0("restrictleft=", nchar(sort(c(Fbarcodes, Rbarcodes))[1])))
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

    if (!dir.exists("logs")) {
      dir.create("logs")
      }

    if (exists("logs/bbdemux.log") | exists("logs/bbdemuxsterr.log")) {
      file.remove(c("logs/bbdemux.log", "logs/bbdemuxsterr.log"))
    }

    result <- system2(command="java",
                      args = args,
                      stdout="logs/stdout.log",
                      stderr="logs/stderr.log",
                      wait=TRUE)
    now <- date()
    cat(paste0("Executed: ", now, "\n"), file="logs/bbdemux.log", append=TRUE)
    file.append("logs/bbdemuxsterr.log", "logs/stderr.log")
    file.remove(c("logs/stdout.log", "logs/stderr.log"))
  }


  if (nsamples > 1) {
    for (i in 1:nsamples) {
      bbtools_seal(
        install = install, fwd = fwd[i], rev = rev[i], Fbarcodes = Fbarcodes, Rbarcodes = Rbarcodes,
        restrictleft = restrictleft, outpath = outpath, kmer = kmer, hdist = hdist, degenerate = degenerate,
        overwrite = overwrite, threads = threads, mem = mem
      )
    }
  } else if (nsamples == 1) {
    bbtools_seal(install, fwd, rev,
      Fbarcodes = Fbarcodes, Rbarcodes = Rbarcodes, restrictleft = restrictleft,
      outpath = outpath, kmer = kmer, hdist = hdist, degenerate = degenerate,
      overwrite = overwrite, threads = threads, mem = mem
    )
  }
}


# Split interleaved reads -------------------------------------------------

#' Split interleaved reads
#'
#' @param install (Required) Install location for bbmap
#' @param files (Required) Vector of locations of interleaved read files to split
#' @param overwrite (Optional) Default TRUE
#' Option to overwrite existing output files.
#'
#' @return
#' @export
#'
#' @examples
bbtools_split <- function(install = NULL, files, overwrite = FALSE) {
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

    reformat_args <- paste(" -cp ", paste0(install, "/current jgi.ReformatReads "), file, out1, out2, overwrite, collapse = " ")
    system2("java", args = reformat_args)
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
