# Propagate taxonomic assignments to species level ------------------------

#' Propagate taxonomy
#'
#' @param tax
#' @param from
#'
#' @return
#' @export
#'
#' @examples
propagate_tax <- function(tax, from = "Family") {
  .Deprecated(new="seqateurs::propagate_tax", package="seqateurs", old="taxreturn::propagate_tax")
  col.prefix <- substr(colnames(tax), 1, 1) # Assumes named Kingdom, ...

  # Highest level to propagate from
  if (from == "Phylum") (start <- 2)
  if (from == "Class") (start <- 3)
  if (from == "Order") (start <- 4)
  if (from == "Family") (start <- 5)
  if (from == "Genus") (start <- 6)
  if (from == "Species") (start <- 7)

  # Propagate
  for (col in seq(start, ncol(tax))) {
    prop <- is.na(tax[, col]) & !is.na(tax[, col - 1])
    newtax <- tax[prop, col - 1]
    needs.prefix <- !grepl("^[A-z]__", newtax)
    newtax[needs.prefix] <- paste(col.prefix[col - 1], newtax[needs.prefix], sep = "__")
    tax[prop, col] <- newtax
  }
  tax
}

#' summarise_fasta
#'
#' @param x The location of a fasta file or gzipped fasta file.
#' @param label optional, Add an extra column with a label
#' @param origin optional, a table with sequence id numbers and their database origins
#'
#' @return
#' @export
#' @import Biostrings
#' @import stringr
#' @import dplyr
#'
#'
#' @examples
summarise_fasta <- function(x, label=NULL, origin=NULL) {
  if(is.null(origin)){
    out <- Biostrings::fasta.index(x) %>%
      dplyr::mutate(taxid = desc %>%
               stringr::str_remove(pattern="(;)(.*?)(?=$)")  %>%
               stringr::str_remove(pattern="(^)(.*?)(?<=\\|)")) %>%
      dplyr::summarise(nseqs = n(),
                nspecies=n_distinct(taxid),
                mean_length = mean(seqlength),
                q0 = quantile(seqlength, probs=0),
                q25 = quantile(seqlength, probs=0.25),
                q50 = quantile(seqlength, probs=0.5), # Q50 is median
                q75 = quantile(seqlength, probs=0.75),
                q100 = quantile(seqlength, probs=1)
      )

  } else if(is.data.frame(origin) | is_tibble(origin)){
    out <- Biostrings::fasta.index(x) %>%
      dplyr::mutate(taxid = desc %>%
               stringr::str_remove(pattern="(;)(.*?)(?=$)")  %>%
               stringr::str_remove(pattern="(^)(.*?)(?<=\\|)")) %>%
      dplyr::mutate(seqid = desc %>%
               stringr::str_remove(pattern="(\\|)(.*?)(?=$)"))  %>%
      dplyr::left_join(origin, by="seqid") %>%
      dplyr::group_by(origin) %>%
      dplyr::summarise(nseqs = n(),
                nspecies=n_distinct(taxid),
                mean_length = mean(seqlength),
                q0 = quantile(seqlength, probs=0),
                q25 = quantile(seqlength, probs=0.25),
                q50 = quantile(seqlength, probs=0.5), # Q50 is median
                q75 = quantile(seqlength, probs=0.75),
                q100 = quantile(seqlength, probs=1)
      )
  }
  if(is.character(label)) {
    out <- out %>%
      dplyr::mutate(label  = label)
  }
  return(out)
}


# Reformat hierarchy ------------------------------------------------------

#' Reformat hierarchy
#'
#' @param x
#' @param quiet
#' @param ranks
#' @param cores The number of  cores to use for get_ott_lineage
#'
#' @return
#' @export
#' @import ape
#' @import stringr
#' @import tidyr
#' @import tibble
#' @import dplyr
#' @import magrittr
#'
#'
#' @examples
reformat_hierarchy <- function(x, db = NULL, quiet = FALSE,
                               ranks = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), cores=1) {
  time <- Sys.time() # get time
  # Convert to DNAbin
  if (!is(x, "DNAbin")) {
    x <- ape::as.DNAbin(x)
  }
  if (missing(db)) {
    stop("A taxonomic database needs to be provided, generate one with get_ncbi_lineage or get_ott_taxonomy")
  }
  if(attr(db, "type") == "ncbi"){
    # Split current names
    seqnames <- names(x) %>%
      stringr::str_remove(";$") %>%
      stringr::str_split_fixed(";", n = Inf) %>%
      tibble::as_tibble() %>%
      tidyr::separate(col = V1, into = c("acc", "tax_id"), sep = "\\|")%>%
      dplyr::select(1:2, tail(names(.), 1)) %>%
      magrittr::set_colnames(c("acc", "tax_id", "species")) %>%
      dplyr::mutate(tax_id = as.numeric(tax_id))

    # Get lineage from taxid
    lineage <- seqnames %>%
      dplyr::left_join (db %>% dplyr::select(-species), by = "tax_id")  %>%
      tidyr::unite(col = Acc, c(acc, tax_id), sep = "|") %>%
      tidyr::unite(col = names, c(!!ranks), sep = ";")
  } else if(attr(db, "type") == "OTT"){
    lineage <- get_ott_lineage(x, db=db, ranks=ranks, cores=cores) %>%
      tidyr::unite(col = names, c(!!ranks), sep = ";")
  }

  names(x) <- paste(lineage$Acc, lineage$names, "", sep = ";") %>%
  stringr::str_remove(";$")
  time <- Sys.time() - time
  if (!quiet) (message(paste0("Reformatted annotations for ", length(x), " sequences in ", format(time, digits = 2))))
  return(x)
}


# Reformat DADA2 Genus ------------------------------------------------------

# make this just run output hierarchy but removing species rank

#' Reformat DADA2 Genus
#'
#' @param x
#' @param quiet
#' @param ranks
#'
#' @return
#' @export
#'
#' @examples
reformat_dada2_gen <- function(x, db = NULL, quiet = FALSE, ranks = NULL, force=FALSE) {
  if (!is.null(ranks)) {
    ranks <- ranks
  } else if (is.null(ranks)) {
    ranks <- c("kingdom", "phylum", "class", "order", "family", "genus")
    message("No ranks supplied, using default ranks for DADA2 assignTaxonomy: kingdom;phylum;class;order;family;genus")
  }
  x <- reformat_hierarchy(x, db = db, quiet = quiet, ranks = ranks)
  return(x)
}


# Reformat DADA2 Species ----------------------------------------------------

#' Reformat DADA2 Species
#'
#' @param x
#' @param quiet
#'
#' @return
#' @export
#' @import ape
#' @import stringr
#' @import tibble
#' @import magrittr
#' @import dplyr
#'
#' @examples
reformat_dada2_spp <- function(x, quiet = FALSE) {
  time <- Sys.time() # get time
  # Convert to DNAbin
  if (!is(x, "DNAbin")) {
    x <- ape::as.DNAbin(x)
  }
  # Split current names
  seqnames <- names(x) %>%
    stringr::str_replace(";$", "") %>%
    stringr::str_split_fixed(";", n = Inf) %>%
    tibble::as_tibble() %>%
    dplyr::select(1, tail(names(.), 1)) %>%
    magrittr::set_colnames(c("acc", "species")) %>%
    dplyr::mutate(species = stringr::str_replace(species, pattern="_", replacement=" "))

  names(x) <- paste(seqnames$acc, seqnames$species, sep = " ")
  time <- Sys.time() - time
  if (!quiet) (message(paste0("Reformatted annotations for ", length(x), " sequences in ", format(time, digits = 2))))
  return(x)
}


# get_lineage -------------------------------------------------------------

#' Get lineage
#'
#' @param x
#' @param db
#'
#' @return
#' @export
#' @import stringr
#' @import tibble
#' @import dplyr
#' @import tidyr
#'
#' @examples
get_lineage <- function(x, db){
  if(missing(db)){ db <- taxreturn::get_ncbi_lineage()}
  cat("Getting taxonomic lineage from taxids\n")
  lineage <- names(x) %>%
    stringr::str_split_fixed(";", n = 2) %>%
    tibble::as_tibble() %>%
    tidyr::separate(col = V1, into = c("acc", "tax_id"), sep = "\\|") %>%
    dplyr::rename(species = V2) %>%
    dplyr::mutate(tax_id = as.numeric(tax_id))  %>%
    dplyr::left_join (db %>% dplyr::select(-species), by = "tax_id")  %>%
    tidyr::unite(col = V1, c(acc, tax_id), sep = "|")  %>%
    dplyr::rename(Acc = V1)
  return(lineage)
}


# train IDTAXA -----------------------------------------------------------

#' Train IDTAXA
#'
#' @param x
#' @param maxGroupSize
#' @param maxIterations
#' @param allowGroupRemoval
#' @param db
#' @param quiet
#' @param force
#' @param get_lineage Reformat the sequences heirarchially before training
#'
#' @return
#' @export
#' @import Biostrings
#' @import ape
#' @import DECIPHER
#' @import stringr
#'
#' @examples
train_idtaxa <- function(x, maxGroupSize=10, maxIterations = 3,  allowGroupRemoval = TRUE,  get_lineage=FALSE, db = NULL, quiet = FALSE) {
  time <- Sys.time() # get time

  #Reformat to complete taxonomic hierarchy
  if(get_lineage & !is.null(db)){
  seqs <- taxreturn::reformat_hierarchy(x, db=db, quiet=FALSE)
  } else if(get_lineage &  is.null(db)){
    stop("If get_lineage is TRUE, a db needs to be provided")
  } else  (seqs <- x)

  #Remove NA's
  if(length(names(seqs)[str_detect(names(seqs), ";NA;")]) > 1){
    remove <- names(seqs)[str_detect(names(seqs), ";NA;")]
    subset <- seqs[!names(seqs) %in% remove]
    message(paste0(length(seqs) - length(subset)," Sequences with NA's in taxonomy removed"))
    seqs <- subset
  }

  #Convert to DNAstringset for DECIPHER
  seqs <- seqs %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet

  # As taxonomies are encoded in the sequence names rather than a separate file, use:
  taxid <- NULL

  #Add Root Rank
  names(seqs) <- names(seqs)  %>% stringr::str_replace(";[;]*", ";Root;")

  # obtain the taxonomic assignments
  groups <- names(seqs) # sequence names

  # assume the taxonomy begins with 'Root;'
  groups <- gsub("(.*)(Root;)", "\\2", groups) # extract the group label
  groupCounts <- table(groups)
  u_groups <- names(groupCounts) # unique groups

  if (!quiet) {message(paste0(length(u_groups), " unique species"))}

  # Pruning training set
  remove <- logical(length(seqs))
  for (i in which(groupCounts > maxGroupSize)) {
    index <- which(groups==u_groups[i])
    keep <- sample(length(index),
                   maxGroupSize)
    remove[index[-keep]] <- TRUE
  }
  if (!quiet) {message(paste0(sum(remove), " sequences pruned from over-represented groups"))}

  # Training the classifier
  probSeqsPrev <- integer() # suspected problem sequences from prior iteration
  for (i in seq_len(maxIterations)) {
    if (!quiet) {message("Training iteration: ", i, "\n", sep="")}
    # train the classifier
    trainingSet <- DECIPHER::LearnTaxa(seqs[!remove], names(seqs)[!remove], taxid)

    # look for problem sequences
    probSeqs <- trainingSet$problemSequences$Index
    if (length(probSeqs)==0) {
      if (!quiet) {message("No problem sequences remaining.\n")}
      break
    } else if (length(probSeqs)==length(probSeqsPrev) &&
               all(probSeqsPrev==probSeqs)) {
      if (!quiet) {message("Iterations converged.\n")}
      break
    }
    if (i==maxIterations)
      break
    probSeqsPrev <- probSeqs

    # remove any problem sequences
    index <- which(!remove)[probSeqs]
    remove[index] <- TRUE # remove all problem sequences
    if (!allowGroupRemoval) {
      # replace any removed groups
      missing <- !(u_groups %in% groups[!remove])
      missing <- u_groups[missing]
      if (length(missing) > 0) {
        index <- index[groups[index] %in% missing]
        remove[index] <- FALSE # don't remove
      }
    }
  }
  if (!quiet) {message(paste0(sum(remove), " sequences removed"))}
  if (!quiet) {message(paste0(length(probSeqs), " problem sequences remaining"))}

  time <- Sys.time() - time
  if (!quiet) (message(paste0("trained IDTAXA on ", length(x), " sequences in ", format(time, digits = 2))))
  return(trainingSet)
}


# taxonomy_to_newick ------------------------------------------------------

#' tax2tree
#'
#' @param x a DNAbin object or an object coercible to DNAbin
#' @param ranks The taxonomic ranks currently assigned to the names
#' @param summarise Select a taxonomic rank to summarise at
#' @param output The output to return, options are:
#' "phylo" to return an ape phylo object
#' "data.tree" to return a data.tree node object
#' "treedf" to return a simplified tree with taxon summaries in a data frame
#' "newick" to return a newick file for visualisation in other software
#'
#' @return a phylo, data.tree, newick, or treedf object
#' @import data.tree
#' @import ape
#' @import stringr
#' @import dplyr
#' @import tidyr
#' @import magrittr
#' @export
#'
#' @examples
tax2tree <- function(x, ranks = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), summarise = "species", output="treedf"){

  # start timer
  time <- Sys.time()
  # Convert to DNAbin
  if (!is(x, "DNAbin")) {
    x <- ape::as.DNAbin(x)
    if (all(is.na(ape::base.freq(x)))) {stop("Error: Object is not coercible to DNAbin \n")}
  }

  if (!output %in% c("phylo", "treedf", "data.tree", "newick")) { stop("Output must be 'phylo', 'data.tree', 'newick' or 'treedf'")}

  # leave an autodetect for ranks?
  ranks <- stringr::str_to_lower(ranks)
  summarise <- stringr::str_to_lower(summarise)
  groupranks <- ranks[1:match(summarise, ranks)]

  #Get taxonomic lineage and convert to tree
  lineage <- names(x) %>%
    stringr::str_replace(pattern=";$", replacement = "") %>%
    stringr::str_split_fixed(";", n=Inf) %>%
    as.data.frame() %>%
    magrittr::set_colnames(c("Acc", ranks )) %>%
    dplyr::select(-Acc) %>%
    dplyr::group_by_at(groupranks) %>%
    dplyr::summarise(sum = dplyr::n()) %>%
    tidyr::unite(col=pathString, !!groupranks, sep="/") %>%
    dplyr::mutate(pathString = paste0("Root/", pathString)) %>%
    data.tree::as.Node(.)

  if (output=="phylo"){
  out <-  ape::read.tree(textConnection(data.tree::ToNewick(lineage, heightAttribute = NULL)))
  } else if (output=="treedf"){
  out <- data.tree::ToDataFrameTree(lineage, "sum")
  } else if (output=="data.tree"){
    out <- lineage
  } else if (output =="newick"){
    out <- data.tree::ToNewick(lineage, heightAttribute = NULL)
  }

  time <- Sys.time() - time
  message(paste0("Generated a taxonomic tree for ", length(x), " Sequences in ", format(time, digits = 2)))

  return(out)
}

# LCA probs ---------------------------------------------------------------

#' Probabilitys of sharing a rank as a function of sequence identity
#'
#' @param x a DNAbin object or an object coercible to DNAbin
#' @param method The distance matrix computation method to use,
#' accepts "mbed" which computes a distance matrix from each sequence to a subset of 'seed' sequences using the method outlined in Blacksheilds et al (2010).
#' This scales well to big datasets, alternatively "kdist" computes the full n * n distance matrix.
#' @param k integer giving the k-mer size used to generate the input matrix for k-means clustering.
#' @param nstart value passed to  nstart passed to kmeans. Higher increases computation time but can improve clustering accuracy considerably.
#' @param ranks The taxonomic ranks currently assigned to the names
#' @param delim The delimiter used between ranks
#'
#' @return
#' @export
#' @import ape
#' @import stringr
#' @import kmer
#' @import dplyr
#'
#'
#' @examples
lca_probs <- function(x, method="mbed",  k=5, nstart = 20, ranks=c("kingdom", "phylum", "class", "order", "family", "genus", "species"), delim=";"){
  # Convert to DNAbin
  if (!is(x, "DNAbin")) {
    x <- ape::as.DNAbin(x)
  }

  # Replace any trailing delimiters
  names(x) <- names(x) %>%
    stringr::str_remove(";$")
  # Count number of remaining delimiters
  ndelim <- stringr::str_count(names(x)[1], ";")

  #check if delims match ranks
  if(any(ndelim <=1)) {
    stop("Error: Needs to be in hierarchial format first, run get_lineage or get_ott_lineage")
  } else if (!ndelim==length(ranks)){
    stop("Error: Number of delimiters does not match number of ranks")
  }

  # Create distance matrix - could probably do this using random samples
  if(method == "kdist"){
    dist <- as.matrix(kmer::kdistance(x, k=k, method="edgar", residues="DNA"))
  } else if(method == "mbed"){
    dist <- as.matrix(kmer::mbed(x, k=k, residues="DNA")[,])
  }

  #convert to pairwise distances
  xy <- t(combn(colnames(dist), 2))
  pw <- data.frame(xy, dist=(100 - round(dist[xy] * 100)))

  # subset to the values in loop

  sim <- sort(unique(pw$dist), decreasing = TRUE)
  simlist <- vector("list", length=length(sim))
  s=1
  for (s in 1:length(sim)) {

    subsets <- pw %>%
      filter(dist==sim[s])

    df1 <- subsets %>%
      tidyr::separate(X1, into=c("Acc",ranks), sep=";") %>%
      dplyr::select(rev(ranks))

    df2 <- subsets %>%
      tidyr::separate(X2, into=c("Acc",ranks), sep=";") %>%
      dplyr::select(rev(ranks))

    #Get all shared ranks
    logidf <- as.data.frame(df1 == df2)

    #Get lowest common rank
    keepvec <- apply(logidf, 1, which.max)
    rows <- seq(logidf[,1])
    #cols <- seq(logidf[1,])
    #unselect <- matrix(ncol=2,c(rep(rows, length(cols)), sort(rep(cols, length(rows)))))
    select <- matrix(ncol=2,c(rows, keepvec))

    logidf[select] <- "KEEP"
    logidf[!logidf=="KEEP"] <- 0
    logidf[logidf=="KEEP"] <- 1
    logidf <- logidf  %>%
      mutate_all(as.numeric) %>%
      colSums() / length(rows)

    simlist[[s]]  <- data.frame(rank=names(logidf), prob=logidf, sim=sim[s])

  }
  names(simlist) <- sim
  out <- dplyr::bind_rows(simlist) %>%
    dplyr::group_by(rank, sim) %>%
    dplyr::summarise(prob = mean(prob))
  return(out)
}



# gid_to_acc --------------------------------------------------------------


#' Convert NCBI gene ids to accession numbers
#'
#' @param ids A character or numeric vector of NCBI gids
#' @param db The origin database of the gids, default 'nuccore'
#' @param chunksize The size of the chunked searches to conduct.
#' Warning, chunk sizes over 300 can be too big for the NCBI servers.
#' @param multithread Whether multithreading should be used
#' @param quiet (Optional) Print text output
#'
#' @return
#' @export
#' @import future
#' @import furrr
#' @import rentrez
#' @import purrr
#'
#'
#' @examples
gid_to_acc <- function(ids, db="nuccore", chunksize=300, multithread=TRUE, progress=FALSE, quiet=FALSE){
  if(!class(ids) %in% c("character", "numeric")){
    stop("input ids must be a character or numeric vector of NCBI gids")
  }
  time <- Sys.time() # get time

  # split search into chunks
  chunks <- split(ids, ceiling(seq_along(ids)/chunksize))
  if(!quiet){message(paste0("Converting ", length(ids), " NCBI gids to accession numbers in ", length(chunks), " chunks"))}

  # setup multithreading
  if(multithread){
    future::plan(future::multiprocess)
  } else if(!multithread){
    future::plan(future::sequential)
  }

  #Main function
  out <- furrr::future_map(chunks, function(x){
    upload <- rentrez::entrez_post(db=db, id=x)
    dl <- rentrez::entrez_fetch(db = database, web_history = upload, rettype = "acc", retmax = chunksize)
    acc <- readLines(textConnection(dl))
    acc <- acc[!acc==""]
    return(acc)
  }, .progress = progress) %>%
    unlist(use.names=FALSE)

  #finished
  time <- Sys.time() - time
  if (!quiet) {message(paste0("Sucessfully converted ", length(out), " of ", length(ids), " NCBI gids to accession numbers in ", format(time, digits = 2)))}
  return(out)
}

# Accessions from fasta ---------------------------------------------------

#' Genbank accessions from fasta
#'
#' @param x A filepath or vector of filepaths to fasta files
#'
#' @return
#' @export
#' @import Biostrings
#' @import dplyr
#'
#' @examples
acc_from_fasta <- function(x) {
  if(!all(file.exists(x))){
    stop("Input must be a filepath or vector of filepaths to fasta files")
  }
  out <- Biostrings::fasta.index(x) %>%
    dplyr::mutate(acc = desc %>%
                    stringr::str_remove(pattern="\\|.*$")) %>%
    dplyr::pull(acc)

  return(out)
}
