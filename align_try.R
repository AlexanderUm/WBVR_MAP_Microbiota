# Finding of species that are close to each other. 
# The idea is to performe aligmnet of each sequence with the sequence in Database
# and find closest mutch 

#1. Set enviorment 

set.seed(9356)

library(msa)


# 0.3 Optimization for parralel processing 
# 0.3.1 Load required libraries 
library(foreach)
library(doParallel)


# 0.3.2 Register claster 
cores=detectCores()
cl <- makeCluster(cores[1]-1) 
registerDoParallel(cl)

# 0.3.3 Remove temp objects 
rm(list = c("cl", "cores"))

#2. Get test data 
#Sequences of taxa from database1 
ps0.I <- readRDS("C:/Rprojects/Microbiota_MAP/output/processing/ps1.1.rds")
ps0.II <- readRDS("C:/Rprojects/Microbiota_MAP_II_20200417/output/dada2/phyloseq0.rds")

taxt.I <- rownames(ps0.I@tax_table@.Data)


otut.II <- ps0.II@otu_table@.Data

sp.II <- names(otut.II[1, ][otut.II[1, ] > 0])


align.msa <- function(seq1_vect, seq2_vect, seq1_poition, seq2_position) {
  
  seq.c <- c(seq1_vect[seq1_poition], seq2_vect[seq2_position])
  
  alg <- msa::msa(seq.c, type = "DNA")
  
  miss <- str_count(msaConsensusSequence(alg), "\\?")
  
  samm.alig <- c(seq.c, miss)
  
  return(samm.alig)
  
  gc()
}

align.msa(taxt.I, sp.II, 36, 2)




fus <- foreach(i=1:14, .packages = c("msa", "stringr")) %dopar% {align.msa(taxt.I, sp.II, i, 2)}
