###################################
# ASV picking using DADA2 pipeline
###################################

# According to http://benjjneb.github.io/dada2/tutorial.html

# Prepare folder for dada2 output: 
dir.create("output/2_dada2")

dir.create("output/2_dada2/diagnostic_plots")

dir.create("output/2_dada2/objects")



# Load libraries 
source("scr/functions/general/load_abs_install_pkg.R")

load_abs_install_pkg(c("dada2", "phyloseq", "tidyverse"))



# Make a list of fastq files in the directory 
path = "data/Dada2_input/"

files <- list.files(path)

files <- files[grep(".fastq", files)]



# Separate forward and revers reads into different variables 
fnFs <- sort(list.files(path, pattern="R1.fastq", full.names = TRUE))

fnRs <- sort(list.files(path, pattern="R2.fastq", full.names = TRUE))


# Make a variable containing only files names 
sample.names <- unique(gsub("_R1.fastq|_R2.fastq", "", files))


# Check reads quality by graphing 
#   Plots are saved in "output/dada2/diagnostic_plots"
forwReads <- plotQualityProfile(fnFs[1:4])

revReads <-plotQualityProfile(fnRs[1:4])

ggsave("output/2_dada2/diagnostic_plots/forward_reads_qualty.pdf", plot = forwReads )

ggsave("output/2_dada2/diagnostic_plots/revers_reads_qualty.pdf", plot = revReads )

# Prepare directory and names for filtered files 
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))

filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names

names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen=c(280, 240), # Truncet bad reads at the end of the reads 
                     # however keep in mind that reads should overlap at least 20bp
                     # in our case amplicon is 444 nucleotide long, so total length should be 464 at least. 
                     trimLeft = c(17, 21), # Trim primers from the begining of sequence                       
                     maxN=0, 
                     maxEE=c(1,1), 
                     truncQ=2, 
                     rm.phix=TRUE,
                     compress=TRUE, 
                     multithread=TRUE) #only for Linux  



# Learn errors in sequences for forward and reverse reads 
errF <- learnErrors(filtFs, multithread=TRUE, nbases = 2e8)

errR <- learnErrors(filtRs, multithread=TRUE, nbases = 2e8)

save(list = c("errF", "errR"), file="output/2_dada2/objects/error_rate_FR.RData")


errFplot <- plotErrors(errF, nominalQ=TRUE)

errRplot <- plotErrors(errR, nominalQ=TRUE)


# Save plots in output/dada2/diagnostic_plots/
ggsave("output/2_dada2/diagnostic_plots/error_rate_F.pdf", plot = errFplot)

ggsave("output/2_dada2/diagnostic_plots/error_rate_R.pdf", plot = errRplot)

# Dereplication - combine identical sequences into unique with corresponding number 
#                    of reads 
derepFs <- derepFastq(filtFs, verbose=FALSE)

derepRs <- derepFastq(filtRs, verbose=FALSE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names

names(derepRs) <- sample.names

# Save dereplicated objects 
save(derepFs, file="output/2_dada2/objects/derepFs.RData")

save(derepRs, file="output/2_dada2/objects/derepRs.RData")



# Samples inference 
#      *addiotion of (..., pool=TRUE) argument could increase sensetivity for rear variates 
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)

dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Back up generated objects 
save(dadaFs, file="output/2_dada2/objects/dadaFs.RData")

save(dadaRs, file="output/2_dada2/objects/dadaRs.RData")


# Merge paired reads 
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Back up generated object 
save(mergers, file="output/2_dada2/objects/mergers.RData")

# Construct sequence table
seqtab <- makeSequenceTable(mergers)

# Check visually reads distribution by their merged length
reads.per.seqlen <- tapply(colSums(seqtab), factor(nchar(getSequences(seqtab))), sum)

df <- data.frame(length=as.numeric(names(reads.per.seqlen)), count=reads.per.seqlen)

readsDisr <- ggplot(data=df, aes(x=length, y=count)) + geom_col() + theme_bw()

# Save the plot 
ggsave("output/2_dada2/diagnostic_plots/leng_reads_distr.pdf", plot = readsDisr)

# Trim reads outside of reads distribution borders 
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(400, 430)]


# Remove chimeras (will be majority of variants)
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)

# Save ASV as R object and csv table 
save(seqtab.nochim, file = "output/2_dada2/objects/seqtab_nochim.RData")

write.csv(seqtab.nochim, file = "output/2_dada2/seqtab_nochim.csv")

# Check how many sequenses get through 
getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), 
               sapply(mergers, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. 
# replace sapply(dadaFs, getN) with getN(dadaFs)

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")

rownames(track) <- sample.names

write.csv(track, file = "output/2_dada2/diagnostic_plots/track.csv")

load("output/2_dada2/objects/seqtab_nochim.RData")

# Taxonomy assignment  
#  Assign taxonomy using classical method 
taxa <- assignTaxonomy(seqtab.nochim, 
                       "resourses/silva_nr_v138_train_set.fa.gz",
                       multithread = TRUE)

# Save taxonomy table as csv and R object
write.csv(x = taxa, "output/2_dada2/taxa.csv")

save(taxa, file ="output/2_dada2/objects/taxa.RData")


# Read in metadata for samples 
metadata_final <- read.csv("output/1_shedding_analysis/metadata_f.csv")

rownames(metadata_final) <- metadata_final$NewId

rownames(seqtab.nochim) <- gsub("\\.gz", "", rownames(seqtab.nochim))


# Combine taxonomic table, ASV count table and metadata into a pholoseq object. 
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), 
               sample_data(metadata_final), 
               tax_table(taxa))

# Save phyloseq object 
save(ps, file = "output/2_dada2/phyloseq0.RData")
saveRDS(ps, "output/2_dada2/phyloseq0.rds")
