
###################################################
# Taxa significantly contributing to classification 
###################################################


# Set enviornment 
#################

set.seed(9356947)

source("scr/functions/general/load_abs_install_pkg.R")

load_abs_install_pkg(c("phyloseq", "tidyverse", "rfPermute", "ComplexHeatmap"))

dir.create("output/7_discriminatory_tax")


# Build random forest model with permutations
# rfPremute package
#############################################

# Prepare data for RF model 
# * from the source
###########################

ps.tf.css.01 <- readRDS("output/3_filtering_mormalization/ps_tf2_css.RDS")

source("scr/functions/data_for_rf.R")

rf.data.css <- data_for_rf(phyloseq = ps.tf.css.01, 
                              class.column = 'WeightedScoreII', 
                              remove.taxa.prev.less.than = 1, 
                              return.df = TRUE)

rf.data.css$WeightedScoreII <- as.numeric(as.character(rf.data.css$WeightedScoreII))

colnames(rf.data.css)[colnames(rf.data.css) %in% "WeightedScoreII"] <- "Shedder"

rf.data.css$Shedder <- as.factor(ifelse(rf.data.css$Shedder < 0.51, "Low", "High"))

# Adjust input data 
# Current row names give an error 

# rf.data.perm <- rf.data.css

# rownames(rf.data.perm) <- paste0("Sample_", 1:nrow(rf.data.css))


# Build the random forest model with optimal parameters
# and find taxa significantly contributing to classification 
############################################################

set.seed(43957)
rf.perm.obj <- rfPermute(Shedder ~ ., 
                         data = rf.data.perm, 
                         importance=TRUE,
                         proximity=TRUE, 
                         mtry = 565,
                         ntree= 15001, 
                         nrep = 999,
                         num.cores = 32)

save(rf.perm.obj, file = "output/plots/7_discriminatory_tax/srf_perm_obj.RData")


# Extract features that are significantly affecting mean accuracy and the Gini index
####################################################################################

pval.rf.perm <- data.frame(rf.perm.obj$pval)

pval.rf.perm.s0 <- pval.rf.perm[pval.rf.perm$MeanDecreaseAccuracy.scaled < 0.05, ]

pval.rf.perm.s1 <- pval.rf.perm.s0[pval.rf.perm.s0$MeanDecreaseGini.scaled < 0.05, ]

# Plot abundance of significantly contributing taxa 
########################################################

# Prepare data 
abund.heat <- rf.data.perm[, colnames(rf.data.perm) %in% rownames(pval.rf.perm.s1)]

# Order rows by status 
abund.heat <- abund.heat[order(rf.data.perm$Shedder), ]

# Order columns by abundance  
abund.heat <- abund.heat[,order(colSums(abund.heat)) ]

# Prepare data for colored bar 
col.stat <- c("steelblue", "gold3")

names(col.stat) <- c("High", "Low")

heat.ant.0 <- data.frame(as.character(rf.data.perm$Shedder))

heat.ant.0$rnam <- rownames(rf.data.perm)

rownames(heat.ant.0) <- rownames(rf.data.perm)

heat.ant.0 <- heat.ant.0[rownames(abund.heat), ]

heat.ant <- data.frame(heat.ant.0[,1])

colnames(heat.ant) <- "Shedding"

rownames(heat.ant) <- heat.ant.0$rnam

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

pdf("output/7_discriminatory_tax/heatmap_sigt_relative.pdf", width = 6, height = 3.5, paper='special')
g.heat.abound
dev.off()

# Box plot - local function
###########################

box_sig <- function(data.matrix, heat.ant, prev) {
    
    #Prepare data 
    if (prev == TRUE) {
        
        High <- (as.vector(colSums(data.matrix[heat.ant$Shedding %in% "High", ]))/
                  length(heat.ant$Shedding[heat.ant$Shedding %in% "High"]))*100

        Low <- (as.vector(colSums(data.matrix[heat.ant$Shedding %in% "Low", ]))/
                  length(heat.ant$Shedding[heat.ant$Shedding %in% "Low"]))*100
        
    } else {
        
        High <- as.vector(colMeans(data.matrix[heat.ant$Shedding %in% "High", ]))
        Low <- as.vector(colMeans(data.matrix[heat.ant$Shedding %in% "Low", ]))   
        
    }


    d <- data.frame(High = High, Low = Low)

    # Add samples names as column names 
    d$Taxa <- colnames(data.matrix)[1:ncol(data.matrix)]

    # Add column with an indication if taxa decrease or increase in the groups  
    d$colr <- ifelse(d$High <= d$Low, "Increase", "Dicrease")

    # Melt the dataframe 
    d.m <- reshape2::melt(d) 

    # Convert column Taxa into factor 
    d.m$Taxa <- factor(d.m$Taxa)

    # Add value for a proper log transformation 
    if (prev == FALSE) {
        d.m$value <- d.m$value + 0.01
    }


# Plot the results 
    if (prev == TRUE) {
        
        p.box <- ggplot(d.m, aes(x = variable, y = value)) +   
                  geom_point() + 
                  geom_line(aes(group = Taxa, color=colr), size = 0.8, alpha=0.75) +   
                  geom_boxplot(alpha=0.01) + theme_bw() + theme(legend.position = "none") + 
                  ylab("Prevalence (%)") + 
                  xlab("Status")
        
    } else {
        
        p.box <- ggplot(d.m, aes(x = variable, y = log(value))) +   
                  geom_point() + 
                  geom_line(aes(group = Taxa, color=colr), size = 0.8, alpha=0.75) +   
                  geom_boxplot(alpha=0.01) + theme_bw() + theme(legend.position = "none") + 
                  ylab("Abundance (Log)") + 
                  xlab("Status")
        
    }
    
return(p.box)
    
}


sig.abund.box <- box_sig(data.matrix = abund.heat, heat.ant = heat.ant, prev = FALSE)

ggsave("output/7_discriminatory_tax/MeanAbandSig.pdf", sig.abund.box, width = 3, height = 5)

ggsave("output/7_discriminatory_tax/MeanAbandSig.png", sig.abund.box, width = 3, height = 5, dpi = 300)


# Plot prevalence difference between groups 
##########################################

# Heatmap data 
prev.heat <- abund.heat

prev.heat[prev.heat > 0] <- 1 

# Plot the heatmap 
set.seed(69734)

g.heat.prev <- ComplexHeatmap::Heatmap(t(prev.heat), 
                                         top_annotation = ha, 
                                         show_column_names = FALSE, 
                                         #km = 2,
                                         show_row_names = FALSE,
                                         #col = my_palette, 
                                         cluster_columns = FALSE,
                                         cluster_rows = TRUE,
                                         name = "Abundance" )

# Save plotted results as figure 
pdf("output/7_discriminatory_tax/heatmap_sig_prev.pdf", width = 6, height = 3.5, paper='special')
g.heat.prev
dev.off()


# Plot boxplot 
sig.prev.box <- box_sig(data.matrix = prev.heat, heat.ant = heat.ant, prev = TRUE)

ggsave("output/7_discriminatory_tax/PrevSig.pdf", sig.prev.box, width = 3, height = 5)

ggsave("output/7_discriminatory_tax/PrevSig.png", sig.prev.box,  width = 3, height = 5, dpi = 400)


# Taxonomy of significantly contributing taxa 
#############################################

# Select only significantly contributing and core taxa 
tax.table <- tax_table(ps.tf.css.01)

sig.tax <- data.frame(tax.table[rownames(tax.table) %in% colnames(abund.heat), ])

# Adjust genus names 
sig.tax.m <- mutate_all(sig.tax, as.character)

rownames(sig.tax.m) <- rownames(sig.tax)

for (i in 1:nrow(sig.tax.m)) {
    
    if (is.na(sig.tax.m$Genus[i])) {
        
        t.na <- sig.tax.m[i,]
        
        t.na <- t.na[!is.na(t.na)]
        
        sig.tax.m$Genus[i] <- paste0("Unknown", "_(", t.na[length(t.na)], ")") 
    }
}

# Summary of information about taxa significant for classification.  

sig.tax.abund <- rf.data.perm[, rownames(sig.tax.m)]

sig.tax.prev <- sig.tax.abund

sig.tax.prev[sig.tax.prev > 0] <- 1 

sig.tax.m$MeanAdunanceLow <- colMeans(sig.tax.abund[rf.data.perm$Shedder %in% "Low", ])

sig.tax.m$MeanAdunanceHigh <- colMeans(sig.tax.abund[rf.data.perm$Shedder %in% "High", ])

sig.tax.m$PrevalenceLow <- colMeans(sig.tax.prev[rf.data.perm$Shedder %in% "Low", ]) * 100

sig.tax.m$PrevalenceHigh <- colMeans(sig.tax.prev[rf.data.perm$Shedder %in% "High", ]) * 100

write.csv(sig.tax.m, "output/7_discriminatory_tax/sig_contr_taxa.csv")


# Identify core of significantly contributing taxa 
###################################################
# A taxa identified as core if present in more than 25% of observed population

# Subset data 

core.prev <- prev.heat[, (colSums(prev.heat)/nrow(prev.heat)) > 0.25]

core.sig.tax <- data.frame(sig.tax.m[colnames(core.prev), ])

core.abund <- abund.heat[, colnames(core.prev)]

write.csv(core.sig.tax, "output/7_discriminatory_tax/core_sig_contr_taxa.csv")


# Preparation of the data for plotting as a heatmap 
core.sig.tax.f <- core.sig.tax[order(as.character(core.sig.tax$Genus)), ]

core.sig.tax.f$Genus <- sub("UCG-005", "Oscillospiraceae_UCG-005", core.sig.tax.f$Genus)

core.sig.tax.f$Genus <- sub("dgA-11_gut_group", "Rikenellaceae_dgA-11", core.sig.tax.f$Genus)

core.sig.tax.f$Genus <- gsub("_gut|_group", "", core.sig.tax.f$Genus)

core.abund.h <- core.abund[, rownames(core.sig.tax.f)]

colnames(core.abund.h) <- make.unique(core.sig.tax.f$Genus)


#PLot and save the heatmap
##########################

core.adund.heat <- ComplexHeatmap::Heatmap(t(core.abund.h), 
                        row_names_side = "left", 
                        top_annotation = ha, 
                        show_column_names = FALSE, 
                        #km = 2,
                        show_row_names = TRUE,
                        #col = my_palette, 
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        name = "Abundance" )

pdf("output/7_discriminatory_tax/heatmap_coresig_relative.pdf", width = 15, height = 8, paper='special')
core.adund.heat
dev.off()


# Box plot - abundance of core significant taxa 
core.sig.abund.box <- box_sig(data.matrix = core.abund.h, heat.ant = heat.ant, prev = FALSE)

ggsave("output/7_discriminatory_tax/AbundCoreSig.pdf", core.sig.abund.box, width = 3, height = 5)

ggsave("output/7_discriminatory_tax/AbundCoreSig.png", core.sig.abund.box,  width = 3, height = 5, dpi = 400)


# Box plot - prevalence of core significant taxa 
core.prev.h <- core.abund.h

core.prev.h[core.prev.h > 0] <- 1

core.sig.prev.box <- box_sig(data.matrix = core.prev.h, heat.ant = heat.ant, prev = TRUE)

ggsave("output/7_discriminatory_tax/PrevCoreSig.pdf", core.sig.prev.box, width = 3, height = 5)

ggsave("output/7_discriminatory_tax/PrevCoreSig.png", core.sig.prev.box,  width = 3, height = 5, dpi = 300)
