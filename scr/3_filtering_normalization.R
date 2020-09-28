
##########################################
# Filtering and normalization of ASV count 
##########################################


# load libraries
set.seed(32426)

source("scr/functions/general/load_abs_install_pkg.R")

load_abs_install_pkg(c("phyloseq", "metagenomeSeq"))


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

#Remove sample that have less than a 1000 observations 
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
