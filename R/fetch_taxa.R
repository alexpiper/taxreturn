#Split the two functions into fetch_GB and fetch_BOLD
#Use a third function - fetch_taxa, that combines these two with the option of selecting either?
#Switch for loops to safely(Purrr) and then a function afterwards ot check if it worked
#Add an option to merge the file or put in seperate files, also add an option to gzip or not

#For the genbank, can we use the raw JSON files to get the taxonomy, rather than requireing a second step with taxonomizr

fetchtaxa <- function(x, database="bold",marker="COI-5P",quiet=FALSE,downto="family",cores=1){
  #Get downstream taxa
  taxon <- bind_rows(downstream(x, db = "ncbi", downto = downto))

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
      cores <- parallel::makeCluster(cores, outfile="")
      junk <- clusterEvalQ(cores, sapply(c("bold","taxize","tidyverse","rentrez"), require, character.only = TRUE)) #Discard result
      para <- TRUE
      stopclustr <- TRUE
    }else{
      para <- FALSE
      stopclustr <- FALSE
    }
  }

  #Define bold search function
  boldSearch <- function(search,marker,quiet){
    time <- Sys.time() # get time
    data <- bold_seqspec(taxon=search, sepfasta=FALSE)
    if(length(data)!=0 & !is.na(data)){
      #Add Domain level to match GenBank download
      data$domain_name <- "Eukaryota"
      data <- subset(data, select=c("processid", "phylum_name","domain_name", "class_name",
                                    "order_name", "family_name", "genus_name",
                                    "species_name","markercode", "nucleotides")) %>%
        na.omit()   %>% #Remove all rows with NA to get rid of insufficiently ID'd specimens
        dplyr::filter(grepl(marker, markercode)) %>% #Remove all sequences for unwanted markers
        dplyr::filter(!grepl("sp.", species_name))

      #Turn into fasta
      bold_seqname <- subset(data, select=c("processid", "domain_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name"))
      bold_seqname <- apply(bold_seqname, 1, paste, collapse=";")
      bold_seqname <- str_replace_all(bold_seqname," ","_")
      data <- as.character(data$nucleotides)

      cat(file=paste0("bold/",search,"_BOLD.fa")) # delete old file
      for (i in 1:length(data)){
        exp <- paste(">", bold_seqname[i], "\n", data[i], "\n", sep="")
        cat(exp, file=paste0("bold/",search,"_BOLD.fa"), append=T)
        gzip(filename=paste0("bold/",search,"_BOLD.fa"))
      }
      time <- Sys.time() - time
      if (quiet==FALSE){
        message(paste0("Downloaded ", length(bold_seqname)," ", search, " Sequences from BOLD ", " in ", format(time, digits=2)))
      }
    } else {
      warning(paste0("No data for ",search," on bold\n"))
    }
    return(NULL)}

  #Define genbank search function - edit this to output the taxid in header
  gbSearch <- function(input, marker,quiet){
    tryCatch(
      {
        searchQ <- paste("(",input, "[ORGN])", " AND (", paste(c(marker), collapse=" OR "), ") AND 1: 2000" ,"[Sequence Length]", sep="")
        message(searchQ)
        search_results <- entrez_search(db   = "nuccore", term = searchQ, retmax=9999999, use_history=TRUE)
        message(search_results$count)


        if(search_results$count!=0 & !is.na(search_results$count)){
          message(paste0(search_results$count," Sequences to be downloaded for: ", input))

          destfile <- paste0("genbank/",input,"_gb.fa")
          cat(file = destfile, sep="") # delete old file

          l <- 1
          start <- 0

          chunks <- length(search_results$ids)/10000
          if (!is.integer(chunks)){chunks <- as.integer(length(search_results$ids)/10000)+1}

          for(l in 1:chunks){

            #dl <- entrez_fetch(db="nuccore", web_history= search_results$web_history, rettype="fasta", retmax=10000, retstart= start)
            #cat(dl, file= destfile, sep=" ", append=T)

            dl <- entrez_fetch(db="nuccore", web_history= search_results$web_history, rettype="fasta", retmax=10000, retstart= start)
            gb <- gbRecord(rcd=textConnection(dl))
            feat <- getFeatures(gb)
            #Use entrez_fetch rettype ="gb"
            #read as biofiles file using gbrecord(rcd=textConnection(dl))

            message("Chunk", l, " of ",chunks, " downloaded\r")
            start <- start + 10000
            Sys.sleep(2.5)


            if (l >= chunks){
              #delete old zipped file
              fn=paste0(destfile,".gz")
              if (file.exists(fn))
                #Delete file if it exists
                file.remove(fn)

              gzip(filename=destfile)
            }
          }
        }}, error=function(e) NULL)
    return(NULL)}


  #Search function
  if(database=="bold"){
    if (!file.exists("bold")){ dir.create("bold")}
    checks <- bold_tax_name(taxon$childtaxa_name)
    search <- checks$taxon[which(!is.na(checks$taxon))]
    search <- if(para){
      parallel::parLapply(cores, search, boldSearch,marker,quiet)
    }else{
      lapply(search, boldSearch,marker,quiet)
    }
    if(para & stopclustr) parallel::stopCluster(cores)
    if(!quiet) message("Done\n")
  }  else if(database=="genbank"){
    if (!file.exists("genbank")){ dir.create("genbank")}
    gbtaxon <- taxon$childtaxa_name
    gbtaxon <- if(para){
      parallel::parLapply(cores, gbtaxon, gbSearch,marker,quiet)
    }else{
      lapply(gbtaxon, gbSearch,marker,quiet)

    }
    if(para & stopclustr) parallel::stopCluster(cores)
    if(!quiet) message("Done\n")

  }
}
