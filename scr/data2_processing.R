
##############################################################
# Step 2 
# DADA2 pipline  - as discribed in tutrial (v 1.8)
# The output file will be a 
##############################################################

# if DADA2 is not istalled it could be doen via folowing piece of code
######################################################################
# 1. Install BiocManager if it is not installed yet
#
#    if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# 2. Install data2 itself 
#
#BiocManager::install("dada2") 
#
######################################################################
#Aditional packages 
#Reshape2
######################################################################

#0. Prepare folder for dada2 output: 
dir.create("output/dada2")
dir.create("output/dada2/diagnostic_plots")
dir.create("output/dada2/objects")

#1. Load libraries 
library("dada2")
library("ggplot2")
library("phyloseq")


#2. Set directory of the project 
setwd("c://Rprojects/Microbiota_MAP_II_20200417")


#3. create variable for the path 
path = "data/processed/raw_sequences/"


#4. Make a list of fastq files in the directory 
files <- list.files(path)
files <- files[grep(".fastq", files)]


#5. Separate forward and revers reads into different variables 
fnFs <- sort(list.files(path, pattern="R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="R2.fastq", full.names = TRUE))


#6. Make variabe contaning only files names 
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


#7. check reads quality by graphing 
#   Plots are saved in "output/dada2/diagnostic_plots"
forwReads <- plotQualityProfile(fnFs[1:4])
revReads <-plotQualityProfile(fnRs[1:4])

ggsave("output/dada2/diagnostic_plots/forward_reads_qualty.pdf", plot = forwReads )
ggsave("output/dada2/diagnostic_plots/revers_reads_qualty.pdf", plot = revReads )

rm(list = c("forwReads", "revReads"))


#8. Prepare directory and names for filtered files 
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


#9. Trim sequnces to remove bad bases 
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
                     multithread=FALSE) # ): Windows 
head(out)


#10. Check reads number before and after filtering by ploting (personal addition, not needed for DADA2)
#10.1 Prepare data for the plot 
out.gg <- as.data.frame(out)
out.gg$SampleID <- rownames(out.gg)
out.gg.m <- reshape2::melt(out.gg)
#10.2 Plot using ggplot package 
filtReads <- ggplot(out.gg.m, aes(x=SampleID, y=value, color=variable, group=SampleID)) + 
             geom_point(size = 3.5) + 
             theme_bw() + 
             geom_line(size = 1) + 
             theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
             ylab("Reads number")
#10.3 Save the plot in utput/dada2/diagnostic_plots/ folder 
ggsave("output/dada2/diagnostic_plots/n_reads_filter1.pdf", plot = filtReads, width = 12 )
#10.4 Clean up 
rm(list = c("filtReads", "out.gg", "out.gg.m"))

 
#11. Learn errors in sequenses for forward and reverse reads 
#11.1 Learning itself 
errF <- learnErrors(filtFs, multithread=TRUE, nbases = 2e8)
errR <- learnErrors(filtRs, multithread=TRUE, nbases = 2e8)
#11.1.1 Save  objects with error rates in output data
save(list = c("errF", "errR"), file="output/dada2/objects/error_rate_FR.RData")

#11.2 Plot estimated error rate 
errFplot <- plotErrors(errF, nominalQ=TRUE)
errRplot <- plotErrors(errR, nominalQ=TRUE)
#11.2.1 Save plots in utput/dada2/diagnostic_plots/
ggsave("output/dada2/diagnostic_plots/error_rate_F.pdf", plot = errFplot)
ggsave("output/dada2/diagnostic_plots/error_rate_R.pdf", plot = errRplot)
#11.3 Clean up
rm(list = c("errFplot", "errRplot"))
gc()


#12. Dereplication - combine identical sequenses into unique with corresponding number 
#                    of reads 
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
#12.1 Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
#12.2 Save dereplicated objects 
save(derepFs, file="output/dada2/objects/derepFs.RData")
save(derepRs, file="output/dada2/objects/derepRs.RData")


#13. Samples inference 
#      *addiotion of (..., pool=TRUE) argument could increase sensetivity for rear variates 
#13.1 inference step 
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#13.2 Back up generated objects 
save(dadaFs, file="output/dada2/objects/dadaFs.RData")
save(dadaRs, file="output/dada2/objects/dadaRs.RData")

#13.3 Check number of AVS and basic quality 
dadaFs[[1]]
dadaRs[[1]]


#14.1 Merge paired reads 
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

#14.2 Back up generated object 
save(mergers, file="output/dada2/objects/mergers.RData")


#15.1 Construct sequence table
seqtab <- makeSequenceTable(mergers)

#15.2 Check wisually reads distibution by their merged length
reads.per.seqlen <- tapply(colSums(seqtab), factor(nchar(getSequences(seqtab))), sum)
df <- data.frame(length=as.numeric(names(reads.per.seqlen)), count=reads.per.seqlen)
readsDisr <- ggplot(data=df, aes(x=length, y=count)) + geom_col() + theme_bw()
#15.2.1 Save the plot 
ggsave("output/dada2/diagnostic_plots/leng_reads_distr.pdf", plot = readsDisr)

#15.3 trim reads outside of reads distribution borders 
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(400,430)]


#16. Remove chemeras (will be majority of variants)
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)

#16.1 Save ASV as R object and csv table 
save(seqtab.nochim, file = "output/dada2/objects/seqtab_nochim.RData")
write.csv(seqtab.nochim, file = "output/dada2/seqtab_nochim.csv")


#17. Check how many sequenses get through 
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
write.csv(track, file = "output/dada2/diagnostic_plots/track.csv")


#18. Taxonomy assigment 
#18.1 Assign taxonomy using classical method 
taxa <- assignTaxonomy(seqtab.nochim, "resource/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
#18.1.1 Save taxonomy table as csv and R object
write.csv(x = taxa, "output/dada2/taxa.csv")
save(taxa, file ="output/dada2/objects/taxa.RData")

#18.2 Assign taxonomy to species level using exact match 
taxa.species <- addSpecies(taxa, "resource/silva_species_assignment_v132.fa.gz")


#19. Read in metadata for samples 
load("output/metadata.R")
metadata_final <- meta.clean.bII

rownames(metadata_final) <- metadata_final$Plate.position




#20. Combine taxonomic table, ASV count table and metadata into a pholoseq object. 
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(metadata_final), 
               tax_table(taxa))

#20.1 Save phyloseq object 
save(ps, file = "output/dada2/phyloseq0.RData")
saveRDS(ps, "output/dada2/phyloseq0.rds")

#Clean up the 
rm(list = ls())
gc()
