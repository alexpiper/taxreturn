
#' boldSearch
#'
#' @param x
#' @param marker
#' @param quiet
#' @param output
#' @param file
#' @param compress
#' @param dir
#'
#' @import bold
#' @import tidyverse
#' @import rentrez
#' @import aphid
#' @import insect
#' @import biofiles
#' @import Biostrings
#' @import ape
#' @import taxizedb
#' @import stringr
#' @import parallel
#'
#' @return
#' @export
#'
#' @examples
boldSearch <- function(x,marker="COI-5P",quiet=FALSE,output="h",file=NULL,compress=FALSE,dir=NULL){

  #function setup
  time <- Sys.time() # get time
  if (!output %in% c("h","binom","gb","bold","gb-binom")){ stop(paste0(output, " has to be one of: 'h','binom','bold', 'gb' or 'gb-binom', see help page for more details"))}
  if (marker=="COI-5P")  {if (!quiet) (cat("Using default marker 'COI-5P' \n"))}
  if(is.null(dir)){dir="bold"}
  if (!file.exists(dir)){ dir.create(dir)}
  if (is.null(file)){
    file=paste0(getwd(),"/",dir,"/",x,"_",marker,"_bold.fa")
    if (!quiet) (message(paste0("No input file given, saving output file to: ",file)))
  }
  if (file.exists(file)) {file.remove(file)}
  if (stringr::str_detect(file,".gz")){compress=TRUE}

  #Bold search
  data <- bold::bold_seqspec(taxon=x, sepfasta=FALSE)
  if(length(data)!=0 && !class(data)=="logical"){
    data <- data %>%
      dplyr::na_if("") %>%
      dplyr::filter(grepl(marker, markercode)) %>% #Remove all sequences for unwanted markers
      dplyr::mutate(domain_name = "Eukaryota") %>%
      dplyr::filter(!is.na(species_name))

    if(nrow(data)!=0){
    #Hierarchial output
    if (output=="h") {
      data <- subset(data, select=c("sampleid", "domain_name",
                                    "phylum_name", "class_name",
                                    "order_name", "family_name",
                                    "genus_name", "species_name",
                                    "nucleotides")) %>%
        na.omit() %>%
        tidyr::unite("name",c("sampleid", "domain_name",
                       "phylum_name", "class_name",
                       "order_name", "family_name",
                       "genus_name", "species_name"),
              sep=";") %>%
        dplyr::mutate(name = name %>%
                 stringr::str_replace_all(pattern=" ",replacement="_") %>%
                 trimws(which="both")
        )

    #Binomial output
    } else if (output=="binom") {
      data <- subset(data, select=c("sampleid", "species_name", "nucleotides")) %>%
        na.omit() %>%
        tidyr::unite("name",c("sampleid", "species_name"), sep=";") %>%
        dplyr::mutate(name = name %>%
                 stringr::str_replace_all(pattern=" ",replacement="_") %>%
                 trimws(which="both")
        )

    #BOLD taxID output
    } else if (output=="bold") {
      data <- subset(data, select=c("sampleid", "species_taxID", "nucleotides")) %>%
        na.omit() %>%
        tidyr::unite("name",c("sampleid","species_taxID"), sep=";") %>%
        dplyr::mutate(name = name %>%
                 stringr::str_replace_all(pattern=" ",replacement="_") %>%
                 trimws(which="both")
                 )

    #Genbank taxID output
    } else if (output=="gb") {
      data <- subset(data, select=c("sampleid", "species_name", "nucleotides")) %>%
        na.omit() %>%
        dplyr::mutate(species_name = trimws(species_name, which="both")) %>%
        dplyr::mutate(gb= taxizedb::name2taxid(data$species_name)) %>%
        tidyr::unite("name",c("sampleid","gb"),sep="|")

    #gb-binom
    } else if (output=="gb-binom") {
    data <- subset(data, select=c("sampleid", "species_name", "nucleotides")) %>%
      na.omit() %>%
      dplyr::mutate(species_name = trimws(species_name, which="both"))

    ids <- taxizedb::name2taxid(data$species_name,out_type="summary")
    #Add exception handling for duplicated taxon names
    if (any(duplicated(ids$name_txt))){

      dupname <- ids$name_txt[ duplicated(ids$name_txt)]
      dup <- ids[ids$name_txt %in% dupname,] %>%
        dplyr::mutate(correct = FALSE)
      class <- taxizedb::classification(dup$tax_id)

      for (i in 1:length(class)){
          dup$correct[i] <- any(stringr::str_detect(class[[i]]$name,pattern="Insecta"))
      }

    filt <- dup$tax_id[which(dup$correct==FALSE)]
    ids <- ids %>%
      dplyr::filter(!tax_id %in% filt) %>%
      dplyr::rename(species_name = name_txt)

    data <- data %>%
      dplyr::left_join(ids,by='species_name') %>%
      dplyr::rename(gb = tax_id) %>%
      tidyr::unite("name",c("sampleid","gb"),sep="|") %>%
      tidyr::unite("name",c("name","species_name"),sep=";")

    } else if (!any(duplicated(ids$name_txt))){
    data <- data %>%
      dplyr::mutate(gb = taxizedb::name2taxid(data$species_name)) %>%
      tidyr::unite("name",c("sampleid","gb"),sep="|") %>%
      tidyr::unite("name",c("name","species_name"),sep=";")
    }
    }

    #Output fASTA
    #iupac <- paste0(paste(names(Biostrings::IUPAC_CODE_MAP),collapse="|"),"|-")

    #Problem -some bold sequences contain an Ionisine
    data <- data %>%
      dplyr::filter(!stringr::str_detect(nucleotides, pattern="I"))

    seqs <- Biostrings::DNAStringSet(data$nucleotides)
    names(seqs) <- data$name
    if (compress==TRUE){Biostrings::writeXStringSet(seqs,file,format="fasta",compress="gzip",width=5000)
    } else if(compress==FALSE){Biostrings::writeXStringSet(seqs,file,format="fasta",width=5000)}

    #Done message
    time <- Sys.time() - time
    if (!quiet) (message(paste0("Downloaded ", length(seqs)," ", x, " Sequences from BOLD ", " in ", format(time, digits=2))))
  } } else {
    warning(paste0("No species level data for ",x," on bold\n"))
  }
  invisible(NULL)
}



# Genbank fetching function -----------------------------------------------

#' Genbank search function
#'
#' @param x
#' @param marker
#' @param quiet
#' @param output
#' @param minlength
#' @param maxlength
#' @param file
#' @param compress
#' @param dir
#'
#'
#' @import bold
#' @import tidyverse
#' @import rentrez
#' @import aphid
#' @import insect
#' @import biofiles
#' @import Biostrings
#' @import ape
#' @import taxizedb
#' @import stringr
#' @import parallel
#'
#' @return
#' @export
#'
#' @examples
gbSearch <- function(x, marker="COI", quiet=FALSE,output="h",minlength=1, maxlength=2000,file=NULL,compress=FALSE,dir=NULL){

  #function setup
  time <- Sys.time() # get time
  if (!output %in% c("h","binom","gb","bold","gb-binom")){ stop(paste0(output, " has to be one of: 'h','binom','bold', 'gb' or 'gb-binom', see help page for more details"))}
  if (marker=="COI")  {if (!quiet) (cat("Using default marker 'COI' \n"))}
  if (is.null(file)){
    if(!is.null(dir)){file=paste0(getwd(),"/",dir,"/",x,"_",marker,"_gb.fa")
    } else (file=paste0(getwd(),"/","genbank/",x,"_",marker,"_gb.fa"))
    if (!quiet) (message(paste0("No input file given, saving output file to: ",file)))
    if (!file.exists("genbank")){ dir.create("genbank")}
  }
  if (file.exists(file)) {file.remove(file)}
  if (stringr::str_detect(file,".gz")){compress=TRUE}

  tryCatch(
    {
      #Genbank Search
      searchQ <- paste("(",x, "[ORGN])", " AND (", paste(c(marker), collapse=" OR "), ") AND ", minlength,":", maxlength,"[Sequence Length]", sep="")
      search_results <- rentrez::entrez_search(db   = "nuccore", term = searchQ, retmax=9999999, use_history=TRUE)

      if(search_results$count!=0 & !is.na(search_results$count)){
        if (!quiet) (message(paste0(search_results$count," Sequences to be downloaded for: ", searchQ)))

        l <- 1
        start <- 0

        chunks <- length(search_results$ids)/10000
        if (!is.integer(chunks)){chunks <- as.integer(length(search_results$ids)/10000)+1}

        for(l in 1:chunks){

          #Fetch gb flat files
          dl <- rentrez::entrez_fetch(db="nuccore", web_history= search_results$web_history, rettype="gb", retmax=10000, retstart= start)
          gb <- biofiles::gbRecord(rcd=textConnection(dl))

          #Hierarchial output
          if (output=="h") {
          lineage <- biofiles::getTaxonomy(gb) %>%
            str_split_fixed(pattern=";",n=Inf) %>%
            trimws(which="both") %>%
            as_tibble()         %>%
            dplyr::mutate(Species = biofiles::getOrganism(gb)) %>%
            dplyr::mutate(Genus = str_replace(V15, pattern="[.]", replacement="")) %>%
            tidyr::unite("names", c("V1","V2","V4","V6","V10","V14","Genus","Species"), sep=";") %>%
            dplyr::mutate(names = str_replace(names,pattern=" ", replacement="_"))
          names <- paste0(names(biofiles::getSequence(gb)),";",lineage$names)

            #Genbank taxID output
          } else if (output=="gb") {
              cat_acctax <- function(x){
                taxid <- purrr::map(x,biofiles::dbxref, "taxon")
                tax_chr <- purrr::map_chr(taxid, function(y){as.character(y[1,1])})
                attributes(tax_chr) <- NULL
                taxout <- paste0(attributes(taxid)$names,"|",tax_chr)
                return(taxout)
              }
          names <- cat_acctax(gb)
          } else if (output=="binom") {

            names<- paste0(names(biofiles::getSequence(gb)),";",biofiles::getOrganism(gb))
          } else if (output=="gb-binom"){
            cat_acctax <- function(x){
              taxid <- purrr::map(x,biofiles::dbxref, "taxon")
              tax_chr <- purrr::map_chr(taxid, function(y){as.character(y[1,1])})
              attributes(tax_chr) <- NULL
              taxout <- paste0(attributes(taxid)$names,"|",tax_chr)
              return(taxout)
            }
            names <- paste0(cat_acctax(gb),";",biofiles::getOrganism(gb))

          }

          #Output FASTA
          seqs <- biofiles::getSequence(gb)
          names(seqs) <- names
          if (compress==TRUE){Biostrings::writeXStringSet(seqs,file,format="fasta",compress="gzip",width=5000)
          } else if(compress==FALSE){Biostrings::writeXStringSet(seqs,file,format="fasta",width=5000)}


          if (!quiet) (message("Chunk", l, " of ",chunks, " downloaded\r"))
          start <- start + 10000
          Sys.sleep(2.5)
          if (l >= chunks){
            time <- Sys.time() - time
            if (!quiet) (message(paste0("Downloaded ", length(seqs)," ", x, " Sequences from Genbank ", " in ", format(time, digits=2))))
          }
        }
      }}, error=function(e) NULL)
  invisible(NULL)
}


# Fetchseqs wrapper function ----------------------------------------------


#' Fetchseqs wrapper function
#'
#' @param x
#' @param database
#' @param marker
#' @param downstream
#' @param downto
#' @param quiet
#' @param output
#' @param minlength
#' @param maxlength
#' @param dir
#' @param compress
#' @param cores
#'
#'
#'@import bold
#' @import tidyverse
#' @import rentrez
#' @import aphid
#' @import insect
#' @import biofiles
#' @import Biostrings
#' @import ape
#' @import taxizedb
#' @import stringr
#' @import parallel
#'
#'
#' @return
#' @export
#'
#' @examples
fetchSeqs <- function(x,database, marker="COI", downstream=FALSE,downto="family", quiet=TRUE, output="h", minlength=1, maxlength=2000,dir=NULL,compress=FALSE, cores=1){
  #Setup parallel
  if(inherits(cores, "cluster")){
    para <- TRUE
    stopclustr <- FALSE
  }else if(cores == 1){
    para <- FALSE
    stopclustr <- FALSE
  }else{
    navailcores <- parallel::detectCores()
    if(identical(cores, "autodetect")) cores <- navailcores - 1
    if(cores > 1){
      # if(cores > navailcores) stop("Number of cores is more than number available")
      if(!quiet) cat("Multithreading with", cores, "cores\n")
      cores <- parallel::makeCluster(cores, outfile="out.txt")
      junk <- parallel::clusterEvalQ(cores, sapply(c("bold","taxizedb","tidyverse","rentrez","Biostrings","biofiles"), require, character.only = TRUE)) #Discard result
      para <- TRUE
      stopclustr <- TRUE
    }else{
      para <- FALSE
      stopclustr <- FALSE
    }
  }

  if (is.null(dir)){
    dir=database
    if (!quiet) (message(paste0("No input dir given, saving output file to: ",dir)))
  }
  if (!file.exists(dir)){ dir.create(dir)}


  if (downstream !=FALSE){
    if(!quiet) cat(paste0("Getting downstream taxa to the level of: ", downto,"\n"))
    taxon <- dplyr::bind_rows(taxizedb::downstream(x, db = "ncbi")) %>%
      dplyr::filter(rank == stringr::str_to_lower(!!downto))

    if (nrow(taxon)>0) {taxon <- taxon$childtaxa_name} else (taxon=x)
    if(!quiet) cat(paste0(length(taxon), " downstream taxa found\n"))
  } else (taxon = x)

  #If taxon is a vector of names
  if(database=="genbank" && length(taxon)>0){
    #Multithread
    taxon <- if(para){
      parallel::clusterExport(cl=cores, varlist=c("taxon", "gbSearch", "quiet", "dir","output","minlength","maxlength","compress"), envir=environment())

      parallel::parLapply(cores, taxon, gbSearch, marker, quiet, file=NULL, output, minlength, maxlength, compress)
    }else{
      #Single core
      lapply(taxon, gbSearch, marker, quiet, file=NULL, output, minlength, maxlength, compress)
    }

    } else if (database=="bold" && length(taxon)>0){
      #Check taxon names exist on bold
      if(!quiet) cat("Checking validity of taxon names for BOLD search\n")
      checks <- bold::bold_tax_name(taxon)
      bold_taxon <- checks$taxon[which(!is.na(checks$taxon))]
      if(!quiet) cat(paste0(length(bold_taxon)," of ", length(taxon), " taxa found to be valid names on BOLD\n"))

      #Multithread
      bold_taxon <- if(para){
        parallel::parLapply(cores, bold_taxon, boldSearch, marker, quiet, file=NULL, output, compress)
      }else{
        #Single core
        lapply(bold_taxon, boldSearch, marker, quiet, file=NULL, output, compress)
      }
    #If taxon is a single name
    }else if (database=="genbank" && length(taxon)==0){
      bold::boldSearch(taxon,marker, quiet, file=NULL, output, compress)
    }else if (database=="bold" && length(taxon)==0){
      bold::boldSearch(taxon,marker, quiet, file=NULL, output, compress)
    }

  #Close clusters
  if(para & stopclustr) parallel::stopCluster(cores)
  if(!quiet) message("Done\n")
invisible(NULL)
}


