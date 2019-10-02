## Get ranked lineage

library(tidyverse)
ranked_lineage <- read_tsv("../rankedlineage.dmp",
                           col_names = c("tax_id","tax_name","species","genus","family","order","class","phylum","kingdom","superkingdom"),
                           col_types=("i-c-c-c-c-c-c-c-c-c-"))


get_ranked_lineage <- function (db = "NCBI", synonyms = FALSE)
    {
      if (!identical(db, "NCBI")) {
        stop("Only the NCBI taxonomy database is available in this version\n")
      }
      tmp <- tempdir()
      fn <- "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
      download.file(fn, destfile = paste0(tmp, "/tmp.tar.gz"))
      message("Extracting data\n")
      test <- untar(tarfile = paste0(tmp, "/tmp.tar.gz"), exdir = tmp)
      if (!identical(test, 0L))
        stop(cat(test))
      message("Building data frame\n")

      #Read data frame
      lin <- read_tsv(paste0(tmp, "/rankedlineage.dmp"),
                      col_names = c("tax_id","tax_name","species","genus","family","order","class","phylum","kingdom","superkingdom"),
                      col_types=("i-c-c-c-c-c-c-c-c-c-"))

      #Remove synonyms
      if (synonyms== FALSE) {
        x <- scan(file = paste0(tmp, "/names.dmp"), what = "", sep = "\n",
                  quiet = TRUE)
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

## data frame to newick

## recursion function
traverse <- function(a,i,innerl){
  if(i < (ncol(df))){
    alevelinner <- as.character(unique(df[which(as.character(df[,i])==a),i+1]))
    desc <- NULL
    if(length(alevelinner) == 1) (newickout <- traverse(alevelinner,i+1,innerl))
    else {
      for(b in alevelinner) desc <- c(desc,traverse(b,i+1,innerl))
      il <- NULL; if(innerl==TRUE) il <- a
      (newickout <- paste("(",paste(desc,collapse=","),")",il,sep=""))
    }
  }
  else { (newickout <- a) }
}

## data.frame to newick function
df2newick <- function(df, innerlabel=FALSE){
  alevel <- as.character(unique(df[,1]))
  newick <- NULL
  for(x in alevel) newick <- c(newick,traverse(x,1,innerlabel))
  (newick <- paste("(",paste(newick,collapse=","),");",sep=""))
}

df <- data.frame(x=c('A','A','B','B','B'), y=c('Ab','Ac','Ba', 'Ba','Bd'), z=c('Abb','Acc','Bad', 'Bae','Bdd'))
myNewick <- df2newick(df, TRUE)

library(ape)
mytree <- read.tree(text=myNewick)
plot(mytree)


