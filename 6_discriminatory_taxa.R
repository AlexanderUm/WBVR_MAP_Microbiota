###################################################
# Taxa significantly contributing to classification 
###################################################


# Set seed 
set.seed(9356947)

# Load required libraries 
source("scr/functions/general/load_abs_install_pkg.R")

load_abs_install_pkg(c("phyloseq", "tidyverse", "rfPermute", "ComplexHeatmap", "gmodels"))

# Create output directory 
dir.create("output/6_discriminatory_taxa")

# 1. Build random forest model with permutations
#    rfPremute package
################################################

# Read in filtered and normalized phyloseq object 
ps.tf.css.01 <- readRDS("output/3_filt_norm_phyl/ps_tf2_css.RDS")

# Load custom function to format data from phyloseq for RF analysis 
source("scr/functions/data_for_rf2.R")

# Format data for RF regression analysis
rf.data.css <- data_for_rf2(phyloseq = ps.tf.css.01, 
                              class.column = 'WeightedScoreII', 
                              remove.taxa.prev.less.than = 1)

# Change name of the "WeightedScoreII" to "Shedder"
colnames(rf.data.css)[colnames(rf.data.css) %in% "WeightedScoreII"] <- "Shedder"

# Define all samples with Shedding score less than 0.51 as "Low" Shedders 
rf.data.css$Shedder <- as.factor(ifelse(rf.data.css$Shedder < 0.51, "Low", "High"))

# Copy "rf.data.css" object into a new object 
rf.data.perm <- rf.data.css

# Change row names 
rownames(rf.data.perm) <- paste0("Sample_", 1:nrow(rf.data.css))


# 2. Build the random forest model with optimal parameters
#    and find taxa significantly contributing to classification 
############################################################

set.seed(43957)
rf.perm.obj <- rfPermute(Shedder ~ ., 
                         data = rf.data.perm, 
                         importance=TRUE,
                         proximity=TRUE, 
                         mtry = 565,
                         ntree= 151, 
                         nrep = 9,
                         num.cores = 12)

# Save the resulted object 
#save(rf.perm.obj, file = "output/6_discriminatory_taxa/srf_perm_obj.RData")


# 3. Extract features that are significantly affecting mean accuracy and the Gini index
#######################################################################################

# Extract P values from RF Premute object in as a dataframe 
pval.rf.perm <- data.frame(rf.perm.obj$pval)

# Subset taxa significantly contributing to MeanDecreaseAccuracy (p < 0.05)
pval.rf.perm.s0 <- pval.rf.perm[pval.rf.perm$MeanDecreaseAccuracy.scaled < 0.05, ]

# Subset taxa significantly contributing to MeanDecreaseGini (p < 0.05)
pval.rf.perm.s1 <- pval.rf.perm.s0[pval.rf.perm.s0$MeanDecreaseGini.scaled < 0.05, ]


# 4. Visualize abundance of significantly contributing taxa as a heatmap 
########################################################################

# Extract data about normalized count of significantly contributing taxa 
abund.heat <- rf.data.perm[, colnames(rf.data.perm) %in% rownames(pval.rf.perm.s1)]

# Order rows by shedding status of corresponding cows  
abund.heat <- abund.heat[order(rf.data.perm$Shedder), ]

# Order columns by total abundance
abund.heat <- abund.heat[,order(colSums(abund.heat)) ]

# Prepare data for colored bar 
col.stat <- c("steelblue", "gold3")

# Adjust names in bar data 
names(col.stat) <- c("High", "Low")

# Extract data about shedding status 
heat.ant.0 <- data.frame(as.character(rf.data.perm$Shedder))

# Add Cow IDs as a column for shedding status 
heat.ant.0$rnam <- rownames(rf.data.perm)

# Add Cow IDs as row names  
rownames(heat.ant.0) <- rownames(rf.data.perm)

heat.ant.0 <- heat.ant.0[rownames(abund.heat), ]

heat.ant <- data.frame(heat.ant.0[,1])

# Adjust column name of heat annotation df 
colnames(heat.ant) <- "Shedding"

# Add Cow IDs as row names for heat annotation df 
rownames(heat.ant) <- heat.ant.0$rnam

# Combine annotation 
ha =  HeatmapAnnotation(df = heat.ant, col = list(Shedding = col.stat))

# Plot the heatmap 
set.seed(69734)
g.heat.abound <- ComplexHeatmap::Heatmap(t(abund.heat), 
                        top_annotation = ha, 
                        show_column_names = FALSE, 
                        #km = 2,
                        show_row_names = FALSE,
                        #col = my_palette, 
                        cluster_columns = FALSE,
                        cluster_rows = TRUE,
                        name = "Abundance" )

# Save the heatmap plot 
pdf("output/6_discriminatory_taxa/Figure_5A.pdf", width = 6, height = 3.5, paper='special')
g.heat.abound
dev.off()


# 5. Create a local boxplot function 
####################################

box_sig <- function(data.matrix, heat.ant, prev) {
    
    # Prepare data 
    # If ploting prevalence plot (prevalence data supplyed)
    if (prev == TRUE) {
        
        # Calculate prevalence in the "High" shedding group by finding sum for each taxa in all samples
        #           and dividing by the number of samples
        High <- (as.vector(colSums(data.matrix[heat.ant$Shedding %in% "High", ]))/
                  length(heat.ant$Shedding[heat.ant$Shedding %in% "High"]))*100
    
        # Calculate prevalence in the "High" shedding group by finding sum for each taxa in all samples
        #           and dividing by the number of samples 
        Low <- (as.vector(colSums(data.matrix[heat.ant$Shedding %in% "Low", ]))/
                  length(heat.ant$Shedding[heat.ant$Shedding %in% "Low"]))*100
     
    # If plotting abundance (abundance data supplied)
    } else {
        
        # Find taxa mean abundance in the "High" shedding group 
        High <- as.vector(colMeans(data.matrix[heat.ant$Shedding %in% "High", ]))
        
        # Find taxa mean abundance in the "Low" shedding group 
        Low <- as.vector(colMeans(data.matrix[heat.ant$Shedding %in% "Low", ]))   
        
    }

    # Create a datatable 
    d <- data.frame(High = High, Low = Low)

    # Add samples names as column names 
    d$Taxa <- colnames(data.matrix)[1:ncol(data.matrix)]

    # Add column with an indication if taxa decrease or increase in the groups  
    d$colr <- ifelse(d$High <= d$Low, "Increase", "Dicrease")

    # Melt the dataframe 
    d.m <- reshape2::melt(d) 

    # Convert column "Taxa" into a factor 
    d.m$Taxa <- factor(d.m$Taxa)

    # Add value for a proper log transformation if abundance is plotted
    if (prev == FALSE) {
        d.m$value <- d.m$value + 0.01
    }


    # Plot the results for prevalence mode 
    if (prev == TRUE) {
        
        p.box <- ggplot(d.m, aes(x = variable, y = value)) +   
                  geom_point() + 
                  geom_line(aes(group = Taxa, color=colr), size = 0.8, alpha=0.75) +   
                  geom_boxplot(alpha=0.01) + theme_bw() + theme(legend.position = "none") + 
                  ylab("Prevalence (%)") + 
                  xlab("Status")
    
    # Plot the results for abundance mode (log transformed)     
    } else {
        
        p.box <- ggplot(d.m, aes(x = variable, y = log(value))) +   
                  geom_point() + 
                  geom_line(aes(group = Taxa, color=colr), size = 0.8, alpha=0.75) +   
                  geom_boxplot(alpha=0.01) + theme_bw() + theme(legend.position = "none") + 
                  ylab("Abundance (Log)") + 
                  xlab("Status")
        
    }
  
# Return the boxplot    
return(p.box)
    
}



# 6. Plot abundance difference of taxa significantly contributes to classification 
#         between "High" and "Low" shedding groups 
#####################################################################################

# Generate abundance boxplot using function above 
sig.abund.box <- box_sig(data.matrix = abund.heat, heat.ant = heat.ant, prev = FALSE)

# Save the plot 
ggsave("output/6_discriminatory_taxa/Figure_5B.pdf", sig.abund.box, width = 3, height = 5)
ggsave("output/6_discriminatory_taxa/Figure_5B.png", sig.abund.box, width = 3, height = 5, dpi = 300)


# 7. Plot prevalence difference of taxa significantly contributes to classification 
#         between "High" and "Low" shedding groups 
####################################################################################

# Copy heatmap data in an object 
prev.heat <- abund.heat

# Convert into prevalence 
prev.heat[prev.heat > 0] <- 1 

# Generate prevalence boxplot using function above 
sig.prev.box <- box_sig(data.matrix = prev.heat, heat.ant = heat.ant, prev = TRUE)

# Save the plot 
ggsave("output/6_discriminatory_taxa/Figure_5C.pdf", sig.prev.box, width = 3, height = 5)
ggsave("output/6_discriminatory_taxa/Figure_5C.png", sig.prev.box,  width = 3, height = 5, dpi = 400)


# 8. Format taxonomy of significantly contributing taxa 
#######################################################

# Extract a taxonomic table from the phyloseq object  
tax.table <- tax_table(ps.tf.css.01)

# Subset significantly contributing taxa 
sig.tax <- data.frame(tax.table[rownames(tax.table) %in% colnames(abund.heat), ])

# Convert every column into character class 
sig.tax.m <- mutate_all(sig.tax, as.character)

# Create vector with prefixes for taxonomic level
phy.ind <- c("k_", "p_", "c_", "o_", "f_", "g_")

# Add prefixes to taxonomic levels 
for (i in 1:ncol(sig.tax.m)) {
    
    sig.tax.m[, i] <- paste0(phy.ind[i], sig.tax.m[, i])
    
}

# Add row names 
rownames(sig.tax.m) <- rownames(sig.tax)

# Add column GenusID 
sig.tax.m$GenusID <- sig.tax.m$Genus

# Format GenusID column 
# "i" is a row number 
for (i in 1:nrow(sig.tax.m)) {
        
        # Subset a row 
        t.na <- sig.tax.m[i,]
        
        # Grep objects that are not "NA"
        t.na <- t.na[!grepl("NA", t.na)]
        
        # Grep objects that are not "-"
        t.na <- t.na[!grepl("-", t.na)]
        
        # Select last object (last identified taxonomic level) and add 
        #         ASV 
        sig.tax.m$GenusID[i] <- paste0(t.na[length(t.na)], "-", "ASV") 
}


# 9. Summary of information about taxa significant for classification.
######################################################################

# Subset data about abundance of significantly contributing taxa 
#        from the data used for RF model with permutations
sig.tax.abund <- rf.data.perm[, rownames(sig.tax.m)]

# Copy to make a prevalence object 
sig.tax.prev <- sig.tax.abund

# Convert abundance into prevalence in the prevalence object
sig.tax.prev[sig.tax.prev > 0] <- 1 

# Calculate mean abundance in the "Low" shedding group 
sig.tax.m$MeanAdunanceLow <- colMeans(sig.tax.abund[rf.data.perm$Shedder %in% "Low", ])

# Calculate mean abundance in the "High" shedding group 
sig.tax.m$MeanAdunanceHigh <- colMeans(sig.tax.abund[rf.data.perm$Shedder %in% "High", ])

# Calculate prevalence expressed in % in the "Low" shedding group
sig.tax.m$PrevalenceLow <- colMeans(sig.tax.prev[rf.data.perm$Shedder %in% "Low", ]) * 100

# Calculate prevalence expressed in % in the "High" shedding group 
sig.tax.m$PrevalenceHigh <- colMeans(sig.tax.prev[rf.data.perm$Shedder %in% "High", ]) * 100

# Make an empty vector for significance testing of differences in 
#      abundance between "Low" and "High" shedding group
wt.prev <- c()

# Test difference significance (loop)
# "i" number of a column 
for (i in 1:ncol(sig.tax.abund)) {
   
   # Perform Wilcoxon test  
   pw <- wilcox.test(sig.tax.abund[rf.data.perm$Shedder %in% "Low", i], 
                     sig.tax.abund[rf.data.perm$Shedder %in% "High", i])
    
   # Combine obtained P values into a vector   
   wt.prev <- c(wt.prev, pw$p.value)  
}

# Add P values to the summary table 
sig.tax.m$WilcoxAbundPval <- wt.prev 

# Add adjusted P values to the summary table 
sig.tax.m$WilcoxAbundQvalFDR <- p.adjust(wt.prev, method = "fdr")

# Add information about contribution to the classification model
rf.imp.tax <- rf.perm.obj$importance

# Subset only significantly contributing taxa 
rf.imp.tax <- rf.imp.tax[rownames(rf.imp.tax) %in% rownames(sig.tax.m), ]

# Combine data from both tables 
sig.tax.f <- cbind(sig.tax.m, rf.imp.tax[, c(3:4)])

# Make GenusID unique
sig.tax.f$GenusID <- make.unique(sig.tax.f$GenusID)

# Save the final table 
write.table(sig.tax.f, "output/6_discriminatory_taxa/sig_contr_taxa.txt")

# 10. Plot Mean Decrease Accuracy and Gini indexes
###################################################

# Prepare data 
# Subset data for plotting from table with information about significant contributing taxa
sig.p.d <- sig.tax.f[, c("Phylum", "GenusID", "MeanDecreaseAccuracy", "MeanDecreaseGini")]

# Order taxa by mean decrease accuracy 
sig.p.d <- sig.p.d[order(sig.p.d$MeanDecreaseAccuracy), ]

# Subset only top 20 contributing taxa 
sig.p.d <- sig.p.d[(nrow(sig.p.d) - 20):nrow(sig.p.d), ]

# Melt (gather) data into long format 
sig.p.dm <- gather(sig.p.d, Index, Value, -Phylum, -GenusID)

# Arrange appropriate order of level in GenusID
sig.p.dm$GenusID <- factor(sig.p.dm$GenusID, levels=sig.p.d$GenusID)

# Plot and save 
mda.plot <- ggplot(sig.p.dm, aes(x = GenusID, y = Value, fill = Phylum)) + 
            geom_bar(stat = "identity", width = 0.6, color = "black") +
            facet_grid(.~ Index, scales = "free") + 
            theme_bw() + 
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
            ylab("") + 
            xlab("AVSs with the corresponding genus names") +
            scale_fill_brewer(palette = "Dark2") + 
            coord_flip()

ggsave(filename = "output/6_discriminatory_taxa/Figure_4.pdf", plot = mda.plot, width = 8, height = 5.5)
ggsave(filename = "output/6_discriminatory_taxa/Figure_4.jpg", plot = mda.plot, width = 8, height = 5.5, dpi = 400)

# 11. Differences in abundance and prevalence between "Low" & "High" shedding groups
####################################################################################
# Will gather summary and statistics used in the manuscript text into a single object (list)

all.summ <- list()

all.summ[["abund.dif"]] <- summary(abs(sig.tax.m$MeanAdunanceLow - sig.tax.f$MeanAdunanceHigh))

all.summ[["prev.dif"]] <- summary(abs(sig.tax.m$PrevalenceLow - sig.tax.f$PrevalenceHigh))

all.summ[["abund.dif.test"]] <- wilcox.test(sig.tax.m$MeanAdunanceLow, sig.tax.f$MeanAdunanceHigh)

all.summ[["prev.dif.test"]] <-wilcox.test(sig.tax.m$PrevalenceLow, sig.tax.f$PrevalenceHigh)

all.summ[["abund.low.ci"]] <- ci(sig.tax.m$MeanAdunanceLow)

all.summ[["abund.high.ci"]] <- ci(sig.tax.f$MeanAdunanceHigh)

