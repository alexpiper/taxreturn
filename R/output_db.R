

# Reformat DADA2 Genus ------------------------------------------------------

#make this just run output heirarchy but removing species rank

reformat_dada2_gen <- function(x,quiet=TRUE,ranks=NULL){
    if(!is.null(ranks)){
      ranks = ranks
    } else if (is.null(ranks)){
      ranks=c("superkingdom","kingdom","class","order","family","genus")
      message("No ranks supplied, using default ranks for DADA2 assignTaxonomy: superkingdom;kingdom;class;order;family;genus")
    }
  x <- reformat_heirarchy(x,quiet=quiet,ranks=ranks)
  return(x)
}
# Reformat DADA2 Species ----------------------------------------------------

reformat_dada2_spp <- function(x,quiet=TRUE){
  time <- Sys.time() # get time
  #Convert to DNAbin
  if(!is(x,"DNAbin")){ x <- ape::as.DNAbin(x)}
  seqnames <- stringr::str_split_fixed(names(x),pattern=";",n=2) %>%
    as_tibble()
  names(x) <- paste(seqnames$V1,seqnames$V2,"",sep=" ")
  return(x)
}



# Reformat heirarchy ------------------------------------------------------

reformat_heirarchy <- function(x,quiet=TRUE,ranks=NULL){
  time <- Sys.time() # get time
  #Convert to DNAbin
  if(!is(x,"DNAbin")){ x <- ape::as.DNAbin(x)}

  if(!is.null(ranks)){
    ranks = ranks
  } else if (is.null(ranks)){
    ranks=c("superkingdom","kingdom","class","order","family","genus","species")
    message("No ranks supplied, using default ranks: superkingdom;kingdom;class;order;family;genus;species")
  }

  #first split names
  seqnames <- names(x) %>%
    stringr::str_split_fixed(";",n=2) %>%
    tibble::as_tibble() %>%
    tidyr::separate(col=V1,into=c("acc","taxid"),sep="\\|") %>%
    dplyr::rename(species = V2) %>%
    mutate(taxid = as.numeric(taxid))

  #Get lineage from taxid
  query <- unique(seqnames$taxid)
  lineage <- insect::get_lineage(query, db=db,simplify=TRUE)

  #Get desired items from lists
  ranklist <- list()
  i=1
  for (i in 1:length(lineage)) {
    line <- as.data.frame(rbind(lineage[[i]]),stringsAsFactors=FALSE)
    if (!is.na(line)[1]){
      line <- line %>% subset(select=which(!duplicated(names(.)))) %>% # drop duplicated `no rank` columns
        dplyr::select(one_of(!!ranks)) %>%   # subset to columns if they exist
        dplyr::mutate(taxid = query[i])              #add query row
      ranklist[[i]] <- line
    } else next
  }

  if ("species" %in% ranks){
    df_lineage <- bind_rows(ranklist) %>%
      rename(species_new = species)
  } else (df_lineage <- bind_rows(ranklist))

  seqnames <- seqnames %>%
    dplyr::left_join(df_lineage,by="taxid") %>%
    tidyr::unite(col=V1,c(acc,taxid),sep="|") %>%
    tidyr::unite(col=V2,c(!!ranks),sep=";")

  names(x) <- paste(seqnames$V1,seqnames$V2,"",sep=";")
  time <- Sys.time() - time
  if (!quiet) (message(paste0("Reformatted ", length(x), " sequences in ", format(time, digits=2))))
  return(x)
}




# Reformat RDP --------------------------------------------------------------


# Output IDTAXA -----------------------------------------------------------

