
<!-- README.md is generated from README.Rmd. Please edit that file -->

# taxreturn

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/alexpiper/taxreturn.svg?branch=master)](https://travis-ci.org/alexpiper/taxreturn)
<!-- badges: end -->

taxreturn is an R package for fetching DNA barcode sequences and
associated taxonomic annotations from public databases such as the
Barcode of Life Database (BOLD) and NCBI GenBank, curating these
sequences and formatting them into training sets compatible with popular
taxonomic classifiers used for metabarcoding and marker gene analysis.

## Installation

This package is still in development and not yet available on CRAN. You
can install development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
#devtools::install_github("alexpiper/taxreturn")
```

# To add to fetchseqs

\*add further filters to the download part - such as binomial only

\#to add to filtseqs resolve synonyms\!

\#add database specific functions Fetching heirarchial structure
Reformat to specific taxonomic database structure rdp format \> DADA2
rdp format \> IDTAXA format \>

``` r
library(taxreturn)
library(Biostrings)
library(tidyverse)

## Fetch sequences from GenBank 
genbank <- fetchSeqs("Trioza", database="genbank",downstream=TRUE,quiet=FALSE, downto="Species", marker="COI OR COI OR COX1 OR COXI", output = "gb-binom",compress=FALSE, cores=3)


## Fetch sequences from BOLD
bold <- fetchSeqs("Trioza", database="bold",downstream=TRUE,quiet=FALSE, downto="Species", marker="COI-5P", output = "gb-binom",compress=FALSE, cores=3)

#read in all fastas and merge
gbSeqs <-  readDNAStringSet(sort(list.files("genbank", pattern = ".fa", full.names = TRUE)))
boldSeqs <-  readDNAStringSet(sort(list.files("bold", pattern = ".fa", full.names = TRUE)))
mergedSeqs <- append(gbSeqs, boldSeqs, after=length(gbSeqs))
uniqSeqs <- mergedSeqs[unique(names(mergedSeqs)),] # Remove those sequnce names that are identical across both databases

#remove non-homologous sequences

#model <- data("model", package="taxreturn")
load("C:/Users/ap0y/Dropbox/R/taxreturn/data/model.rda")
filtered <- clean_seqs(uniqSeqs, model,minscore = 100, cores=2, shave=TRUE,maxNs = 0)

#Here - resolve synonyms and check for those that have changed
#for those that have been changed - retrieve a new NCBI taxID using name2rank

#Resole synonyms

filtered <- resolve_synonyms(filtered,subspecies=FALSE,quiet=FALSE,treat_taxid="ignore",higherrank=FALSE,fuzzy=TRUE)

#Save old names into attributes
attributes(filtered)$oldnames <- names(filtered)
#Get names in format for insect::purge
names(filtered) <- names(filtered) %>%
  str_split_fixed(";",n=2) %>%
  as_tibble() %>%
  pull("V1") 

#filter using insect::purge - Could wrap this in a function for ease of use?
db <- insect::taxonomy(db = "NCBI", synonyms = TRUE)

#get unique names only
filtered <- insect::subset.DNAbin(filtered, subset = !duplicated(names(filtered)))

purged  <- insect::purge(filtered, db = db, level = "species", confidence = 0.8,
                  threshold = 0.99, method = "farthest")

#Restore old names
names(purged) <- attributes(purged)$oldnames

#Prune group sizes down to 5 
pruned <- prune_groups(purged,maxGroupSize = 5, discardby="length",dedup=TRUE, quiet = FALSE)

#remove small sequences
#pruned <- pruned[lengths(pruned)>500]

#Change to RDP taxonomy format
test <- reformat_heirarchy(pruned, quiet=FALSE)

test2 <- reformat_dada2_spp(pruned)


#Trim to primer region using virtualPCR from insect package
amplicon <- virtualPCR(filtseqs, up = "ACWGGWTGRACWGTNTAYCC",down= "ARYATDGTRATDGCHCCDGC",cores=3, rcdown = TRUE, trimprimers = TRUE)
writeFASTA(amplicon,"gb_trimmed.fa")
```

``` r

#build PHMM from midori longest - sequences need to be same length
midori <-  Biostrings::readDNAStringSet("MIDORI_LONGEST_20180221_COI.fasta")
insecta_midori <- as.DNAbin(midori[str_detect(names(midori),pattern=";Insecta;"),])
folmer <- insect::virtualPCR(insecta_midori, up = "TITCIACIAAYCAYAARGAYATTGG",down= "TAIACYTCIGGRTGICCRAARAAYCA",cores=2, rcdown = TRUE, trimprimers = TRUE)
filt <- folmer[lengths(folmer)==658]

#Filtered was then aligned in MAFFT - mafft folmer_insecta_fullength.fa > folmer_insecta_fullength_aligned.fa
#alignment was then manually curated in geneious primer
folmer_curated <-  ape::read.dna("folmer_insecta_fullength_aligned_curated.fa",format="fasta")
model <- aphid::derivePHMM(folmer_curated)
```
