##Could turn this into a function with different inputs!
#
#filter_taxa <- function(file,minlength,maxlength,unique=TRUE,binomials=TRUE, removeterms){
#
#
#}
#
#
#library(insect)
###filter bold
#bold <- read.fasta(file="bold/merged_bold.fa", strip.desc = FALSE, as.string = TRUE)
#bold_names <- getAnnot(bold)
#
###Filter failed taxonomy - NA;NA
#preNA <- length(bold)
#bold <- bold[!grepl("NA;NA|;;|;$ ", bold_names)] #Remove NA's introduced by taxonomizr step, as well as sequences with multiple missing ranks and sequences with final missing rank
#bold_filtered <- cbind(preNA,length(bold))
#
###Filter sequences of innapropriate length
#bold <- bold[which(getLength(bold) >200)] ##Remove all sequences below 200bp
#bold_filtered$smallrem <- length(bold)
#bold <- bold[which(getLength(bold)<2000)] #Remove all sequences above 3000bp
#bold_filtered$bigrem <- length(bold)
#
###Filter duplicate sequences
#bold<- unique(bold)
#bold_filtered$unique <- length(bold)
#
###Filter any further erroneous or insufficiently identified sequences
#
#bold_names <- as_tibble(unlist(getAnnot(bold)))
#bold_names <- bold_names %>%
#  dplyr::filter(!str_detect(value, fixed("sp."))) %>%
#  dplyr::filter(!str_detect(value, fixed("aff."))) %>%
#  dplyr::filter(!str_detect(value, fixed("nr."))) %>%
#  dplyr::filter(!str_detect(value, fixed("cf."))) %>%
#  dplyr::filter(!str_detect(value, fixed("nom."))) %>%
#  dplyr::filter(!str_detect(value, fixed("nud."))) %>%
#  dplyr::filter(!str_detect(value, fixed("environment"))) %>%
#  dplyr::filter(!str_detect(value, fixed("undescribed"))) %>%
#  dplyr::filter(!str_detect(value, fixed("unverified"))) %>%
#  dplyr::filter(!str_detect(value, fixed("uncultured"))) %>%
#  dplyr::filter(!str_detect(value, fixed("unidentif"))) %>%
#  dplyr::filter(!str_detect(value, fixed("Bacterium"))) %>%
#  dplyr::filter(!str_detect(value, fixed("wolbachia"))) %>%
#  dplyr::filter(!str_detect(value, fixed("symbiont"))) %>%
#  dplyr::filter(!str_detect(value, fixed("Bacterium"))) %>%
#  dplyr::filter(!str_detect(value, fixed("NA"))) %>%
#  dplyr::filter(!str_detect(value, fixed("error")))
#
#rm_keywords <- bold_names$value
#
#bold <- bold[getAnnot(bold) %in% rm_keywords]
#bold_filtered$keywords <- length(bold)
#
##Write out filtered fasta
#bold_names <- getName(bold)
#write.fasta(bold, bold_names, paste0("bold_tempfilt1.fa"), as.string=FALSE, nbchar=100)
#
##test
#bold <- readFASTA("bold_tempfilt1.fa")
#
##trim to amplicon - BF1 - BR1
##this is probably better to replace with a PHMM directly for the sequencec
#amplicon <- virtualPCR(bold, up = "ACWGGWTGRACWGTNTAYCC",down= "ARYATDGTRATDGCHCCDGC",cores=3, rcdown = TRUE, trimprimers = TRUE)
#writeFASTA(amplicon,"bold_trimmed.fa")
#
#bold_filtered$trimmed <- length(amplicon)
#
#rm(amplicon)
#rm(bold)
#
