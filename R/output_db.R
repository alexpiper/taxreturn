##-Wrap this all in a function for output_db_type
##Make RDP classifier sets
## specify the path to your file of training sequences:
#seqs_path <- "merged_cleaned.fa"
## read the sequences into memory
#seqs <- readDNAStringSet(seqs_path)
#
## Make species level classifier
#names <- names(seqs) %>%
#  str_replace("_"," ") %>%
#  str_split_fixed(";",n= Inf)
#
#spp_names <- paste0(names[,1]," ",names[,8])
#
#gen_names <- as.tibble(names) %>%
#  subset(select=c(V2,V3,V4,V5,V6,V7))
#gen_names <- apply(gen_names, 1, paste, collapse=";")
#
#
##write out species fasta
#dir.create("reference")
#names(seqs) <- spp_names
#writeXStringSet(seqs,"reference/merged_rdp_species.fa.gz", append=FALSE,
#                compress=TRUE, format="fasta",width=1000)
#
##write out genus fasta
#names(seqs) <- gen_names
#writeXStringSet(seqs,"reference/merged_rdp_genus.fa.gz", append=FALSE,
#                compress=TRUE, format="fasta",width=1000)
#
#
