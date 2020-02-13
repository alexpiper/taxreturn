#' fastqc_install
#'
#' @param url (Optional) Default will search for the latest version
#' URL to retrieve fastQC from.
#' @param dest.dir (Optional)  Default "bin"
#' Directory to install fastQC within.
#' @force Whether existing installs should be forcefully overwritten
#'
#' @return
#' @export
#'
#' @examples
fastqc_install <- function(url, dest.dir = "bin", force = FALSE) {
  if (missing(url)) {

    # find the latest version of fastq
    download_page <- xml2::read_html("http://www.bioinformatics.babraham.ac.uk/projects/download.html")
    link_hrefs <- download_page %>%
      rvest::html_nodes("a") %>%
      rvest::html_attr("href")
    fastqc_href <- grep("fastqc/fastqc.*.zip",
                        link_hrefs, perl = TRUE) %>%
      link_hrefs[.] %>% .[1]
    url <- paste0("http://www.bioinformatics.babraham.ac.uk/projects/",
                  fastqc_href)

  }

  if (!dir.exists(dest.dir)) {
    dir.create(dest.dir) # Create first directory
  }

  # Check if dir exists
  if (dir.exists(paste0(dest.dir, "/fastQC")) && force == FALSE) {
    message("Skipped as FASTQC already exists in directory, to overwrite set force to TRUE")
    return(NULL)
  } else  if (dir.exists(paste0(dest.dir, "/fastQC")) && force == TRUE) {
    unlink(paste0(dest.dir, "/fastQC"), recursive = TRUE) # Remove old version
  }


  destfile <- file.path(dest.dir, basename(url))
  if (exists(destfile)) {
    file.remove(destfile) # Remove old zip file
  }
  httr::GET(url, httr::write_disk(destfile, overwrite=TRUE))

  #unzip file
  utils::unzip(destfile, exdir = dest.dir)
  #Remove download
  file.remove(destfile)
}

#' Run FastQC Tool
#' @description Run FastQC Tool on windows or linux
#' @param fq.dir (Required) Default is the current working directory.
#'    path to the directory containing fastq files.
#' @param qc.dir (Optional) path to the FastQC result directory. If NULL, a directory
#'   named FASTQC is created in the current working directory.
#' @param threads (Optional) the number of threads to be used. Default is 1.
#' @param fastqc.path (Required) path to fastqc program
#' @return A new directory containing the fastqc reports for each sample
#'
#' @import processx
#' @examples
#' \dontrun{
#' # Run FastQC: generates a QC directory
#' fastqc(fq.dir)
#' }
#' @export
fastqc <- function(fq.dir,   qc.dir = NULL, threads = 1, fastqc.path = "bin/FastQC/fastqc")
{
  if(is.null(qc.dir)) qc.dir <- file.path(fq.dir, "FASTQC")
   dir.create(file.path(qc.dir))

   if (.Platform$OS.type == "unix") {
     cmd <- paste0(fastqc.path, " ", fq.dir, "/*  --threads ", threads,  " --outdir ", qc.dir)
     result <- system(cmd)

   } else{
  .fq.dir <- paste0(fq.dir, "/*")
  .threads <- paste0("--threads ", threads)
  .qc.dir <- paste0("--outdir ", qc.dir)

  result <- processx::run(command="perl",
                          args =  c(fastqc.path, .fq.dir, .threads, .qc.dir),
                          echo=TRUE,
                          echo_cmd	= TRUE,
                          spinner=TRUE,
                          windows_verbatim_args=TRUE,
                          error_on_status = FALSE,
                          cleanup_tree = TRUE)

}
  return(result)
}
