## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval=FALSE,
  warning = FALSE,
  error = FALSE,
  message = FALSE
)

## ----setup, echo=FALSE--------------------------------------------------------
#  library(taxreturn)
#  library(Biostrings)
#  library(ape)
#  library(insect)
#  library(tidyverse)

## ----eval=FALSE---------------------------------------------------------------
#  devtools::install_github("alexpiper/taxreturn")
#  library(taxreturn)
#  
#  # Load other necessary packages
#  library(Biostrings)
#  library(ape)
#  library(insect)
#  library(tidyverse)

## ----Download all sequences, eval=FALSE, include=TRUE-------------------------
#  ## Fetch sequences from GenBank by searching for a taxon name
#  fetchSeqs("Scaptodrosophila", database="genbank", out_dir="genbank", downstream=FALSE, marker="COI[GENE] OR COX1[GENE] OR COXI[GENE]", output = "gb-binom", compress=TRUE, force=TRUE, multithread =FALSE)
#  
#  ## OR fetch sequences from a species list
#  spp_list <- readLines("species_list.txt")
#  fetchSeqs(spp_list, database="genbank", out_dir="genbank", quiet=FALSE, output = "gb-binom", marker="COI[GENE] OR COX1[GENE] OR COXI[GENE]", compress=TRUE, force=TRUE, multithread = FALSE)

## ----Download BOLD Sequences, eval=FALSE, include=TRUE------------------------
#  ## Fetch sequences from BOLD by searching for a taxon name
#  fetchSeqs("Scaptodrosophila", database="bold", out_dir="bold", downstream=FALSE, marker="COI-5P", output = "gb-binom", compress=TRUE, force=TRUE, multithread = FALSE)
#  
#  ## OR fetch sequences from a species list
#  spp_list <- readLines("species_list.txt")
#  fetchSeqs(spp_list, database="bold", out_dir="bold", marker="COI-5P", output = "gb-binom", compress=TRUE, force=TRUE, multithread = FALSE)

## ----Download mitochondria, eval=FALSE, include=TRUE--------------------------
#  # Fetch mitochondrial genomes from genbank by searching for a taxon name
#  fetchSeqs("Scaptodrosophila", database="genbank", out_dir="genbank", quiet=FALSE, marker="mitochondria", output = "gb-binom", force=TRUE, compress=TRUE, multithread = FALSE)
#  
#  ## OR fetch sequences from a species list
#  spp_list <- readLines("species_list.txt")
#  fetchSeqs(spp_list, database="genbank", out_dir="genbank", quiet=FALSE, marker="mitochondria", output = "gb-binom", force=TRUE, compress=TRUE, multithread = FALSE)

## ----load model---------------------------------------------------------------
#  #This loads the model into the workspace
#  data("model", package="taxreturn")
#  
#  #See what it looks like
#  print(model)

## ----build PHMM, eval=FALSE---------------------------------------------------
#  # Read in sequence dataset to be used in model training
#  seqs <-  readDNAStringSet("MIDORI_LONGEST_20180221_COI.fasta")
#  
#  # Trim the sequences to the amplified region using a virtual PCR
#  amplicon <- insect::virtualPCR(seqs,
#                                 up = "TITCIACIAAYCAYAARGAYATTGG",  #Forward primer
#                                 down= "TAIACYTCIGGRTGICCRAARAAYCA", #Reverse primer
#                               cores=1, rcdown = TRUE, trimprimers = TRUE)
#  
#  #Only retain amplicons of the appropriate length (in this case 658bp)
#  amplicon_filtered <- amplicon[lengths(amplicon) == 658]
#  model <- aphid::derivePHMM(amplicon_filtered)

## ----clean seqs---------------------------------------------------------------
#  #read in all fastas and merge
#  seqs <- c(readDNAStringSet(list.files("genbank", pattern = ".fa",
#                                        full.names = TRUE)),
#                  readDNAStringSet(list.files("bold", pattern = ".fa",
#                                              full.names = TRUE))
#                  )
#  
#  # Remove those sequnce accessions that are identical across both databases
#  # Using a regex that splits to just the accessions
#  uniqSeqs <- seqs[!duplicated(str_extract(names(seqs), "^.*\\|" )),]
#  
#  #remove non-homologous sequences
#  filtered <- map_to_model(uniqSeqs, model, min_score = 100, min_length = 50, max_N=Inf, max_gap=Inf,
#                           shave=TRUE, check_frame=TRUE, kmer_threshold = 0.3, k=5,
#                           multithread = FALSE, quiet=FALSE)
#  

## ----get mixed clusters-------------------------------------------------------
#  # Load the NCBI taxonomy
#  db <- get_ncbi_taxonomy()
#  
#  #Get db
#  mixed_clusters <- get_mixed_clusters(
#      x = filtered, db=db,
#      rank = "genus",
#      threshold = 0.97,
#      return = "consensus",
#      confidence=0.6, quiet = FALSE)
#  
#  # Get accession numbers to remove
#  rem <- mixed_clusters$Acc
#  
#  purged  <- subset.DNAbin(filtered, subset = !str_replace(names(filtered), "(?:.(?!;))+$", "") %in% rem)

## ----Stop codons--------------------------------------------------------------
#  #SGC4 is invertebrate mitochondrial genetic code
#  codon_filt <- codon_filter(purged, genetic_code = "SGC4", tryrc = TRUE, resolve_draws = "majority")

## ----resolve synoynms---------------------------------------------------------
#  resolved <- resolve_synonyms_ncbi(purged)
#  
#  #Check for differences in names
#  names(resolved)[which(!names(resolved) %in% names(purged))]

## ----prune groups-------------------------------------------------------------
#  #Prune group sizes down to 5, removing all identical sequences first
#  pruned <- prune_groups(resolved, max_group_size = 5, discardby = "length", dedup=TRUE, quiet = FALSE)

## ----trim to primer regions---------------------------------------------------
#  #Trim to primer region using virtualPCR from insect package
#  amplicon <- insect::virtualPCR(pruned, up = "ACWGGWTGRACWGTNTAYCC", down = "ARYATDGTRATDGCHCCDGC", cores = 1, rcdown = TRUE, trimprimers = TRUE)
#  
#  #Check the lengths of the trimmed sequences
#  table(lengths(amplicon))

## ----reformat and output------------------------------------------------------
#  #Load NCBI taxonomy
#  db <- get_ncbi_taxonomy()
#  
#  #Reformat to complete taxonomic heirarchy
#  hierarchy <- reformat_hierarchy(amplicon, db=db, quiet=FALSE)
#  
#  #See if this worked
#  head(names(hierarchy))
#  
#  # Save a zipped fasta file of curated reference database
#  insect::writeFASTA(hierarchy, file = "COI_reference_hierarchial.fa.gz", compress = TRUE)

## ----summary unique-----------------------------------------------------------
#  names(hierarchy) %>%
#    str_split_fixed(";", n=8) %>%
#    as_tibble() %>%
#    tidyr::separate(V1, into=c("Sequences", "taxids")) %>%
#    magrittr::set_colnames(c("Sequences", "tax_ids", "kingdom", "phylum", "class", "order", "family", "genus", "species")) %>%
#    summarise_all(n_distinct)

## ----newick tree--------------------------------------------------------------
#  tree <- tax2tree(hierarchy, output="phylo")
#  write.tree(tree, "tree.nwk")

## ----RDP----------------------------------------------------------------------
#  #Load NCBI taxonomy
#  db <- get_ncbi_taxonomy()
#  
#  #Reformat to Kingdom to genus heirarchy suitable for assigntaxonomy classifier in DADA2
#  dada2_gen <- reformat_dada2_gen(amplicon, db = db, quiet = FALSE)
#  insect::writeFASTA(dada2_gen, file = "COI_reference_dada2gen.fa.gz", compress = TRUE)
#  
#  #Reformat to genus species binomials as suitable for assignSpecies in DADA2
#  dada2_spp <- reformat_dada2_spp(amplicon)
#  insect::writeFASTA(dada2_spp, file = "COI_reference_dada2spp.fa.gz", compress = TRUE)

## ----IDTAXA-------------------------------------------------------------------
#  #Load NCBI taxonomy
#  db <- get_ncbi_taxonomy()
#  
#  # Train IDTAXA
#  training_set <- train_idtaxa(amplicon, max_group_size = 10, max_iterations = 3,  allow_group_removal = TRUE,  get_lineage = TRUE, db = db, quiet = FALSE)
#  
#  #Write out training set
#  saveRDS(training_set, file="trained_idtaxa.rds")

