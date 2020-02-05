# Get_ranked_lineage ------------------------------------------------------

#' Download NCBI taxdump
#'
#' @param db
#' @param synonyms
#'
#' @return
#' @export
#'
#' @examples
get_ranked_lineage <- function(db = "NCBI", synonyms = TRUE, force=FALSE) {
  if (!identical(db, "NCBI")) {
    stop("Only the NCBI taxonomy database is available in this version\n")
  }
  tmp <- tempdir()
  if (force == TRUE | !file.exists(paste0(tmp, "/rankedlineage.dmp")) | !file.exists(paste0(tmp, "/tmp.tar.gz"))) {
    message("Downloading NCBI taxonomy database")
    fn <- "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
    download.file(fn, destfile = paste0(tmp, "/tmp.tar.gz"))
    message("Extracting data\n")
    test <- untar(tarfile = paste0(tmp, "/tmp.tar.gz"), exdir = tmp)
  if (!identical(test, 0L)) {
    stop(cat(test))
  }
  file.remove(paste0(tmp, "/tmp.tar.gz"))
  }
  message("Building data frame\n")

  # Read data frame
  lin <- readr::read_tsv(paste0(tmp, "/rankedlineage.dmp"),
    col_names = c("tax_id", "tax_name", "species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom"),
    col_types = ("i-c-c-c-c-c-c-c-c-c-")
  )

  # Remove synonyms
  if (synonyms == FALSE) {
    x <- scan(
      file = paste0(tmp, "/names.dmp"), what = "", sep = "\n",
      quiet = TRUE
    )
    syn <- x[grepl("synonym", x)]
    syn <- strsplit(syn, split = "\t")
    syn <- sapply(syn, function(s) s[c(1, 3)])
    syn <- as.data.frame(t(syn), stringsAsFactors = FALSE)
    syn[[1]] <- as.integer(syn[[1]])
    colnames(syn) <- c("taxID", "name")
    lin <- lin %>%
      filter(!tax_id %in% syn[[1]])
  }

  message("Done\n")
  return(lin)
}



# ncbi_taxid --------------------------------------------------------------

#' Get ncbi taxid's for a taxon name
#'
#' @param x
#' @param db
#'
#' @return
#' @export
#'
#' @examples
ncbi_taxid <- function(x, db=NULL) {

  if (is.null(db)) { db <- get_ranked_lineage()}
  out <-  as.data.frame(x) %>%
    magrittr::set_colnames("tax_name") %>%
    dplyr::left_join (db, by = "tax_name") %>%
    dplyr::pull(tax_id)

  return(out)
  }


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


# summarise fasta ---------------------------------------------------------


#' summarise_fasta
#'
#' @param x The location of a fasta file or gzipped fasta file.
#' @param label optional, Add an extra column with a label
#' @param origin optional, a table with sequence id numbers and their database origins
#'
#' @return
#' @export
#'
#' @examples
summarise_fasta <- function(x, label=NULL, origin=NULL) {
  if(is.null(origin)){
    out <- fasta.index(x) %>%
      mutate(taxid = desc %>%
               str_replace(pattern="(;)(.*?)(?=$)", replacement="")  %>%
               str_replace(pattern="(^)(.*?)(?<=\\|)", replacement="")) %>%
      summarise(nseqs = n(),
                nspecies=n_distinct(taxid),
                mean_length = mean(seqlength),
                q0 = quantile(seqlength, probs=0),
                q25 = quantile(seqlength, probs=0.25),
                q50 = quantile(seqlength, probs=0.5), # Q50 is median
                q75 = quantile(seqlength, probs=0.75),
                q100 = quantile(seqlength, probs=1)
      )

  } else if(is.data.frame(origin) | is_tibble(origin)){
    out <- fasta.index(x) %>%
      mutate(taxid = desc %>%
               str_replace(pattern="(;)(.*?)(?=$)", replacement="")  %>%
               str_replace(pattern="(^)(.*?)(?<=\\|)", replacement="")) %>%
      mutate(seqid = desc %>%
               str_replace(pattern="(\\|)(.*?)(?=$)", replacement=""))  %>%
      left_join(origin, by="seqid") %>%
      group_by(origin) %>%
      summarise(nseqs = n(),
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
      mutate(label  = label)
  }
  return(out)
}



# Reformat DADA2 Genus ------------------------------------------------------

# make this just run output heirarchy but removing species rank

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
  x <- reformat_heirarchy(x, db = db, quiet = quiet, ranks = ranks)
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
#'
#' @examples
reformat_dada2_spp <- function(x, quiet = FALSE) {
  time <- Sys.time() # get time
  # Convert to DNAbin
  if (!is(x, "DNAbin")) {
    x <- ape::as.DNAbin(x)
  }
  seqnames <- stringr::str_split_fixed(names(x), pattern = ";", n = 2) %>%
    as_tibble()
  names(x) <- paste(seqnames$V1, seqnames$V2, sep = " ")
  time <- Sys.time() - time
  if (!quiet) (message(paste0("Reformatted annotations for ", length(x), " sequences in ", format(time, digits = 2))))
  return(x)
}



# Reformat heirarchy ------------------------------------------------------

#' Reformat heirarchy
#'
#' @param x
#' @param quiet
#' @param ranks
#'
#' @return
#' @export
#'
#' @examples
reformat_heirarchy <- function(x, db = NULL, quiet = FALSE, ranks = NULL, sppsep = "_", force=FALSE) {
  time <- Sys.time() # get time
  # Convert to DNAbin
  if (!is(x, "DNAbin")) {
    x <- ape::as.DNAbin(x)
  }
  if (is.null(db)) {
    message("No taxonomy database provided, downloading from NCBI")
    db <- get_ranked_lineage(force=force)
  }

  if (!is.null(ranks)) {
    ranks <- ranks
  } else if (is.null(ranks)) {
    ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
    message("No ranks supplied, using default ranks: kingdom;phylum;class;order;family;genus;species")
  }

  # Split current names
  seqnames <- names(x) %>%
    stringr::str_split_fixed(";", n = 2) %>%
    tibble::as_tibble() %>%
    tidyr::separate(col = V1, into = c("acc", "tax_id"), sep = "\\|") %>%
    dplyr::rename(species = V2) %>%
    dplyr::mutate(tax_id = as.numeric(tax_id))

  # Get lineage from taxid
  lineage <- seqnames %>%
    dplyr::left_join (db %>% dplyr::select(-species), by = "tax_id")  %>%
    tidyr::unite(col = V1, c(acc, tax_id), sep = "|") %>%
    tidyr::unite(col = V2, c(!!ranks), sep = ";")

  names(x) <- paste(lineage$V1, lineage$V2, "", sep = ";")
  time <- Sys.time() - time
  if (!quiet) (message(paste0("Reformatted annotations for ", length(x), " sequences in ", format(time, digits = 2))))
  return(x)
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
#'
#' @return
#' @export
#'
#' @examples
train_idtaxa <- function(x, maxGroupSize=10, maxIterations = 3,  allowGroupRemoval = TRUE,  db = NULL, quiet = FALSE, force=FALSE) {
  time <- Sys.time() # get time

  #Reformat to complete taxonomic heirarchy
  seqs <- taxreturn::reformat_heirarchy(x, db=db, quiet=FALSE, force=force)

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
    trainingSet <- DECIPHER::LearnTaxa(seqs[!remove],
                             names(seqs)[!remove],
                             taxid)

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
#' @param sim The sequence similarities to cluster at
#' @param k integer giving the k-mer size used to generate the input matrix for k-means clustering.
#' @param nstart value passed to  nstart passed to kmeans. Higher increases computation time but can improve clustering accuracy considerably.
#' @param ranks The taxonomic ranks currently assigned to the names
#' @param delim The delimiter used between ranks
#'
#' @return
#' @export
#'
#' @examples
lca_probs <- function(x, sim=seq(0.9,1,0.01), k=5, nstart = 20, ranks=c("kingdom", "phylum", "class", "order", "family", "genus", "species"), delim=";"){
  # Convert to DNAbin
  if (!is(x, "DNAbin")) {
    x <- ape::as.DNAbin(x)
  }

  simlist <- vector("list", length=length(sim))
  for (s in 1:length(sim)) {
    otus <- kmer::otu(x, k=k, threshold=sim[s] ) %>%
      enframe()  %>%
      tidyr::separate(name, into=c("Acc", ranks, "rep"), sep=";")

    ranklist <- vector("list", length=length(ranks))
    names(ranklist) <- ranks
    for (i in 1:length(ranks)){

      ranklist[[i]] <- otus %>%
        select(c(ranks[i], "value")) %>%
        group_by_all( ) %>%
        dplyr::count() %>%
        group_by(value) %>%             # now required with changes to dplyr::count()
        mutate(prop = prop.table(n)) %>%
        ungroup %>%
        summarise(prop = mean(prop))
    }

    simlist[[s]] <- bind_cols(ranklist) %>%
      set_colnames(ranks) %>%
      mutate(sim = sim[s])
  }
  out <- bind_rows(simlist)
  return(out)
}



