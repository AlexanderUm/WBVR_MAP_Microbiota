
##########################################
# Filtering and normalization of ASV count 
##########################################


# load libraries
set.seed(32426)

source("scr/functions/general/load_abs_install_pkg.R")

load_abs_install_pkg(c("phyloseq", "metagenomeSeq", "tidyverse", "RColorBrewer", "gmodels"))


# Read in phyloseq from DADA2
ps.0 <- readRDS("output/2_dada2/phyloseq0.rds")


# Remove taxa that have count or/end prevalence less than 3 
###########################################################

# Have count less than 3 
ps.0f <- prune_taxa(taxa = colSums(ps.0@otu_table) > 3 , x= ps.0) 

# Have privalence less that 3 
otu_prev <- otu_table(ps.0f)

otu_prev[otu_prev > 1] <- 1

ps_tf1 <- prune_taxa(taxa = colSums(otu_prev) > 3 , x= ps.0f)

# Remove sample that have less than a 1000 observations 
######################################################

o.tab.ps_tf1 <- otu_table(ps_tf1)

ps_tf2 <- prune_samples(samples = rownames(o.tab.ps_tf1)[rowSums(o.tab.ps_tf1) > 1000], ps_tf1)


# Normalize ASVs abundance using CSS as implemented in metagenomSeq 
###################################################################

# Prepare data from phyloseq package 
mg.ps_tf2 <- phyloseq_to_metagenomeSeq(ps_tf2)

# Calculate cumulative statistics 
p <- metagenomeSeq::cumNormStatFast(mg.ps_tf2)

# Normalize count 
mg.c.ps_tf2 <- metagenomeSeq::cumNorm(mg.ps_tf2, p = p)

# Convert metagenomSeq object with normolized count a into dataframe  
css.otu.all <- data.frame(otu_table(metagenomeSeq::MRcounts(mg.c.ps_tf2, 
                        norm = TRUE, log = TRUE), taxa_are_rows = FALSE))



# Adjust names and table format of the otu table 
css.otu.all.m <- as.matrix(t(css.otu.all))

rname <- gsub("[.]", "-", rownames(css.otu.all.m))

rownames(css.otu.all.m) <- gsub("X", "", rname)



# Make a phyloseq object with normolized taxa count and save it 
###############################################################
ps_tf2_css <- ps_tf2

ps_tf2_css@otu_table@.Data <- css.otu.all.m 

dir.create("output/3_filtering_mormalization")
saveRDS(ps_tf2_css, "output/3_filtering_mormalization/ps_tf2_css.RDS")

##################################################
# Extract general information about sequencing
##################################################

# Number of samples per animal 
median(table(ps_tf2_css@sam_data$CowN))

write_csv(x = as.data.frame(table(ps_tf2_css@sam_data$CowN)),
          path = "output/3_filtering_mormalization/samples_number.csv")


# Reads information 
###################
# Number of reads per sample
otus.tab <- otu_table(ps_tf2)

# Total reads 
sum(otus.tab)

# Median, min, max reads 
median(rowSums(otus.tab))

min(rowSums(otus.tab))

max(rowSums(otus.tab))

# Plot reads distribution 
#########################
# Prepare data 
reads.pl.df<- as.data.frame(rowSums(otus.tab))

colnames(reads.pl.df) <- "Reads"

reads.pl.df$SampleID <- rownames(reads.pl.df)

reads.pl.df <- reads.pl.df[order(reads.pl.df$Reads), ]

reads.pl.df$SampleID <- factor(reads.pl.df$SampleID, levels=unique(reads.pl.df$SampleID))


# Plot 
reads.pers.p <- ggplot(reads.pl.df, aes(x = SampleID, y = Reads)) + 
  geom_bar(stat="identity") + 
  theme_bw() + 
  geom_hline(yintercept = mean(reads.pl.df$Reads), color="blue") +
  theme(axis.text.x = element_blank()) + 
  ylab("Reads number") + 
  xlab("Samples") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Save plot
ggsave(filename = "output/3_filtering_mormalization/reads_pesr_sample.pdf", 
       plot = reads.pers.p, width = 10, height = 6)

ggsave(filename = "output/3_filtering_mormalization/reads_pesr_sample.jpg", 
       plot = reads.pers.p, width = 10, height = 6, dpi = 400)


# Phylogenetic composition on phylum level 
##########################################
# Prepare data for plotting 
# Glom phyloseq to Phylum level 
ps.tf2.plot <- tax_glom(ps_tf2, taxrank = "Phylum")

# Transform count to relative abundance 
ps.tf2.plot <- transform_sample_counts(ps.tf2.plot, function(x) x / sum(x) * 100)

# Extract otu table 
phy.plot.d <- data.frame(t(otu_table(ps.tf2.plot)))

# Calculate row sums 
phy.plot.rs <- rowSums(phy.plot.d)

# Add column with Phylum names 
phy.plot.d$Phylum <- as.character(tax_table(ps.tf2.plot)[,"Phylum"])
                                       
# Melt df to long format 
phy.plot.dm <- gather(phy.plot.d, SampleID, Abundance, -Phylum)
                                       
# Order factors in column Phylum by abundance 
phy.plot.dm$Phylum <- factor(phy.plot.dm$Phylum, 
                             levels = c(phy.plot.dm$Phylum[order(phy.plot.rs, decreasing = TRUE)]))


# Prepare colors 
plot.col <- c(brewer.pal(n = 8, name = "Dark2"), replicate('#666666', n = 7))


# Plot as a barplot 
phy.plot <- ggplot(phy.plot.dm, aes(x = SampleID, y = Abundance, fill = Phylum)) + 
        geom_bar(stat = "identity") + 
        scale_fill_manual(values = plot.col) + 
        theme_bw() + 
        theme(axis.text.x = element_blank(), 
              panel.grid = element_blank(), 
              axis.ticks.x = element_blank(), )

# Save plots 
ggsave(filename = "output/3_filtering_mormalization/phylum_plot.pdf", width = 10, height = 5.5)

ggsave(filename = "output/3_filtering_mormalization/phylum_plot.jpg", width = 10, height = 5.5, dpi = 400)


# Abundance summary table 
##########################
# Prepare data 
rownames(phy.plot.d) <- phy.plot.d$Phylum

phy.plot.d2 <- phy.plot.d[, ! colnames(phy.plot.d) %in% "Phylum"]


# Abundance summary and CI
sum.t <- as.data.frame(summary(t(phy.plot.d2)))

sum.t2<- as.data.frame(cbind(as.character(sum.t$Var2), 
              str_split(sum.t$Freq, pattern = ":", simplify = TRUE)))

sum.t3 <- spread(sum.t2, V2, V3)

summ.save <- cbind(sum.t3,  round(t(apply(phy.plot.d2, 1, ci)), 2))

write_csv(x = summ.save, path = "output/3_filtering_mormalization/phylum_abundance.csv")
