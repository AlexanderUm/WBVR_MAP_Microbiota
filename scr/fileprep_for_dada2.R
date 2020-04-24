
#####################################################################
#
#This script will preporecess sequence files provided by BaseClear
#for further analysis using DADA2 pipline                    
#
#####################################################################
# Required pakages (if not installed): 
# stringr, tidyr, R.utils
#####################################################################

#1. set directory of the R project
setwd("C://Rprojects/Microbiota_MAP_II_20200417/")

#2. Copy raw data to data/temp folder 
data.path <- paste0("data/raw/", list.files("data/raw/"))

file.copy(from = data.path, to = "data/temp/")

#3. Unzip raw files in data/temp folder 
unzip(zipfile = data.path, exdir = "data/temp")

#4. Copy sequnce files to data/processed
file.copy(from = "data/temp/raw_sequences", to = "data/processed/", recursive=TRUE)

#5. Clean up data/temp folder 
unlink("data/temp/*", recursive = TRUE)

#6. Compare md5 of files 
#6.1 read in name name of the files 
namesFiles <- list.files(path = "data/processed/raw_sequences/")
#6.2 select only names of fastq.gz files 
namesFilesGz <- namesFiles[grepl(".gz", namesFiles)]

#6.3 Calculate md5 sums 
newfilesMD5 <- c()
# following loop will use 
for (i in namesFilesGz) {
  f.name <- paste("data/processed/raw_sequences", i, sep="/")
  md <- tools::md5sum(f.name)
  names(md) <- i 
  newfilesMD5 <- c(newfilesMD5, md) 
}
rm (list=c("i", "md"))

newfilesMD5df <- data.frame(as.character(newfilesMD5))

newfilesMD5df$sample <- names(newfilesMD5)


### Check if md5 sums correspondes to each other 
# Read in provided by BaseClear md5 files 
md5.path <- paste0("data/processed/raw_sequences/", 
                   namesFiles[grep("md5", namesFiles)])

mdprov <- read.delim(md5.path, header = FALSE)

# Combine them togather 
mdprov <- tidyr::separate(mdprov, col = "V1", into = c("md", "sample"), sep = "  ")


mergedmd5 <- merge(mdprov, newfilesMD5df, by="sample")
# Compare md5 
mergedmd5$md == mergedmd5$as.character.newfilesMD5.


# 7. Rename files
######################

  # following loop will shorten up names of each file 
for (i in namesFilesGz) {
  short_names <- paste(
    stringr::str_extract(i, "[[:upper:]]\\d\\d"), 
    stringr::str_extract(i, "R\\d"), sep="_" 
  )
  short_names <- paste0(short_names, ".fastq.gz")
  short_names <- paste("data/processed/raw_sequences", short_names, sep="/")
  path.names <- paste("data/processed/raw_sequences", i, sep="/")
  file.rename(path.names, short_names)
}

rm (list=c("i", "short_names", "path.names"))

# For some reason utar didn't the comprassion method
####################################################

# List and adjust names of files 
# un.name <- list.files(path = "data/processed/raw_sequences")
# un.name <- paste("data/processed/raw_sequences", 
#                  un.name[grepl(".gz", un.name)], 
#                  sep="/")
# #unzip gz files 
# for (i in un.name) {
#   untar(as.character(i))
# }

# The files are ready for processing using DADA2 pipline 
#
#* will clean up enviorment and memory
unlink("data/processed/raw_sequences/*.gz")
rm(list = ls())
gc()
