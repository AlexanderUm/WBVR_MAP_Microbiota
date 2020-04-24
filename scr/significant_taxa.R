
#
#######################################################
# Identification of significantly contributing to accuracy of 
# random forest model ASV 
#######################################################


# 0 Prepare enviorment 
########################

# 0.1 Set working directory 
setwd("C://Rprojects/Microbiota_MAP/")

# 0.2 Load rfPermute library 
library(rfPermute)
library(ggplot2)
library(ComplexHeatmap)
library(phyloseq)
library(tidyverse)
library(ggtree)
library(tidytree) 
library("phangorn")

# 0.3 Load input file (ASV prepared for RF analysis)
load(file="output/processing/rfdata.Rdata")


# 1. Set seed and run RF with 999 permutations 
set.seed(43957)
rf.perm.obj <- rfPermute(Status ~ ., 
                         data = rfdata, 
                         importance=TRUE,
                         proximity=TRUE, 
                         mtry = 352,
                         ntree=15001, 
                         nrep = 999,
                         num.cores = 8)

# 1.2. Save the resulted object 
#save(rf.perm.obj, file = "output/rf_objects/srf_perm_obj.RData")


# 2. Extract featurs that are significantly affecting mean accuracy and Gini
#############################################################################
pval.rf.perm <- data.frame(rf.perm.obj$pval)
pval.rf.perm.s0 <- pval.rf.perm[pval.rf.perm$MeanDecreaseAccuracy.scaled < 0.05, ]
pval.rf.perm.s1 <- pval.rf.perm.s0[pval.rf.perm.s0$MeanDecreaseGini.scaled < 0.05, ]

rm(list = c("pval.rf.perm", "pval.rf.perm.s0"))

save(pval.rf.perm.s1, file = "output/processing/pval_rf_perm.Rdata")

# 3. Plot abundance of significantly contributing taxa 
########################################################

# 3.1 In a form of heatmap 

# 3.1.1 Prepare data
# 3.1.1.1 Substract from abundance table not significant taxa 
abund.heat <- rfdata[, colnames(rfdata) %in% rownames(pval.rf.perm.s1)]
# 3.1.1.2 Order by status 
abund.heat <- abund.heat[order(rfdata$Status), ]
# 3.1.1.3 Order by total abundance in all samples 
abund.heat <- abund.heat[,order(colSums(abund.heat)) ]

# 3.1.2. Prepare data for colored bar 
col.stat <- c("steelblue", "gold")
names(col.stat) <- c("High", "Low")
heat.ant <- data.frame(as.character(rfdata$Status[order(rfdata$Status)]))
colnames(heat.ant) <- "Status"
rownames(heat.ant) <- rownames(rfdata)
ha =  HeatmapAnnotation(df = heat.ant, col = list(Status= col.stat))

# 3.1.3. Set seed and plot the heatmap 
set.seed(69734)
g.heat.abound <- ComplexHeatmap::Heatmap(t(abund.heat*100), 
                        top_annotation = ha, 
                        show_column_names = FALSE, 
                        #km = 2,
                        show_row_names = FALSE,
                        #col = my_palette, 
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        name = "Abundance" )

pdf("output/plots/heatmap_sigt_relative.pdf", width = 6, height = 4, paper='special')
g.heat.abound
dev.off()


# 3.2 In a form of boxplot 
# 3.2.1 Prepre data 
# 3.2.1.1 Sammarize mean abundance per samples in high and low shedding groups 
High <- as.vector(colMeans(abund.heat[heat.ant$Status %in% "High", ]))*100
Low <- as.vector(colMeans(abund.heat[heat.ant$Status %in% "Low", ]))*100
# 3.2.1.2 Combine into the dataframe 
d <- data.frame(High = High, Low = Low)
# 3.2.1.3 Add samples names as column names 
d$Taxa <- colnames(abund.heat)[1:ncol(abund.heat)]
# 3.2.1.4 Add coulmn with indiation if taxa decrese or increase in the groups 
d$colr <- ifelse(d$High <= d$Low, "Increase", "Dicrease")
# 3.2.1.4 Melt the dataframe 
d.m <- reshape2::melt(d) 
# 3.2.1.4 Convert column Taxa into factor 
d.m$Taxa <- factor(d.m$Taxa)

# 3.2.2 Plot the results 
ggplot(d.m, aes(x = variable, y = value)) +   
  geom_point() + 
  geom_line(aes(group = Taxa, color=colr), size = 0.8, alpha=0.75) +   
  geom_boxplot(alpha=0.01) + theme_bw() + theme(legend.position = "none")

ggsave("output/plots/MeanAbandSig.pdf", width = 3, height = 5)
ggsave("output/plots/MeanAbandSig.png", width = 3, height = 5, dpi = 300)


# 4. Plot prevalence of sinigcantly contributing taxa 
# 4.1 Plot as a heatmap 
# 4.1.1 Prepare presents absence data 
prev.heat <- abund.heat

prev.heat[prev.heat > 0] <- 1 

# 4.1.2 Plot the heatmap 
set.seed(69734)

g.heat.prev <- ComplexHeatmap::Heatmap(t(prev.heat), 
                                         top_annotation = ha, 
                                         show_column_names = FALSE, 
                                         #km = 2,
                                         show_row_names = FALSE,
                                         #col = my_palette, 
                                         cluster_columns = FALSE,
                                         cluster_rows = FALSE,
                                         name = "Abundance" )

# 4.1.3 Save ploted results as figure 
pdf("output/plots/heatmap_sig_prev.pdf", width = 6, height = 4, paper='special')
g.heat.prev
dev.off()


# 4.2 Plot as boxplot
# 4.2.1 Prepre data 
# 4.2.1.1 Sammarize mean abundance per samples in high and low shedding groups 
High <- (as.vector(colSums(prev.heat[heat.ant$Status %in% "High", ]))/
          length(heat.ant$Status[heat.ant$Status %in% "High"]))*100

Low <- (as.vector(colSums(prev.heat[heat.ant$Status %in% "Low", ]))/
          length(heat.ant$Status[heat.ant$Status %in% "Low"]))*100

# 4.2.1.2 Combine into the dataframe 
d <- data.frame(High = High, Low = Low)
# 4.2.1.3 Add samples names as column names 
d$Taxa <- colnames(prev.heat)[1:ncol(prev.heat)]
# 4.2.1.4 Add coulmn with indiation if taxa decrese or increase in the groups 
d$colr <- ifelse(d$High <= d$Low, "Increase", "Dicrease")
# 4.2.1.4 Melt the dataframe 
d.m <- reshape2::melt(d) 
# 4.2.1.4 Convert column Taxa into factor 
d.m$Taxa <- factor(d.m$Taxa)

# 4.2.2 Plot the results 
ggplot(d.m, aes(x = variable, y = value)) +   
  geom_point() + 
  geom_line(aes(group = Taxa, color=colr), size = 0.8, alpha=0.75) +   
  geom_boxplot(alpha=0.01) + theme_bw() + theme(legend.position = "none") + 
  ylab("Prevalence") + 
  xlab("Status")

ggsave("output/plots/PrevSig.pdf", width = 3, height = 5)
ggsave("output/plots/PrevSig.png", width = 3, height = 5, dpi = 300)

# 4.3 Remove objects that were generated during the step 
rm(list = c("abund.heat", "d", "d.m", "g.heat.abound", "g.heat.prev", 
            "ha", "heat.ant", "prev.heat", "col.stat", "High", "Low" ))



#5. Taxonomic assigment of taxa important for classificaiton 
############################################################

#5.1 Prepara data 
#5.1.1 Load phyloseq object (ps1.1)
ps1.1 <- readRDS("output/processing/ps1.1.rds")

#5.1.2 Make a subest from tax_table of sinigicanly contributing taxa 
phy.sig.t <- data.frame(ps1.1@tax_table)[colnames(rfdata) %in% rownames(pval.rf.perm.s1), ]

#5.1.3 Add informaiton about mean decrease accuracy 
phy.sig.t$mda.scal <- pval.rf.perm.s1$MeanDecreaseAccuracy.scaled

#5.1.4 Add information about classificaiton 
phy.sig.t$na.genus <- ifelse(is.na(phy.sig.t$Genus)== "TRUE", "Unassigned Genus", "Assigned Genus")
phy.sig.t$na.genus[is.na(phy.sig.t$Family)] <- "Unassigned Family"

#5.2 Plot results 
#5.2.1 Assigment per phylum 
ggplot(phy.sig.t, aes(x=Phylum, y=mda.scal, fill=na.genus)) + 
  geom_bar(stat = "identity", color="black") + 
  #facet_grid( ~ Phylum, scales = "free", space = "free") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_brewer(palette = "Dark2") + 
  guides(fill=guide_legend(title="ASV taxonomic assigment")) + 
  ylab("Mean Decreas accuracy")

ggsave("output/plots/Sig_f_phylum.pdf", width = 4.5, height = 7)
ggsave("output/plots/Sig_f_phylum.png", width = 4.5, height = 7, dpi = 300)


#5.2.2 Contribution per genus 
ggplot(phy.sig.t, aes(x=Genus, y=mda.scal, fill=Phylum)) + 
  geom_bar(stat = "identity", color="black") + 
  #facet_grid( ~ Family, scales = "free", space = "free") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_brewer(palette = "Set1") + 
  ylab("Mean Decreas accuracy")

ggsave("output/plots/Sig_f_genus.pdf", width = 8.5, height = 4)
ggsave("output/plots/Sig_f_genus.png", width = 8.5, height = 4, dpi = 300)


#5.2.3 Contribution per family 
ggplot(phy.sig.t, aes(x=Family, y=mda.scal, fill=Phylum)) + 
  geom_bar(stat = "identity", color="black") + 
  #facet_grid( ~ Phylum, scales = "free", space = "free") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_brewer(palette = "Set1") + 
  ylab("Mean Decreas accuracy")

ggsave("output/plots/Sig_f_family.pdf", width = 9, height = 4)
ggsave("output/plots/Sig_f_family.png", width = 9, height = 4, dpi = 300)




#6. Prevalence and abundance of identified taxa. 
################################################

#6.1 Normalize the abundance using CSS using metagenomSeq package 
#6.1.1 Prepare data from phyloseq package 
mg.ps1.1 <- phyloseq_to_metagenomeSeq(ps1.1)

#6.1.2 Calculate cumulative statistics 
p <- metagenomeSeq::cumNormStatFast(mg.ps1.1)

#6.1.3 Normalize count 
mg.ps1.1 <- metagenomeSeq::cumNorm(mg.ps1.1, p = p)

#6.1.5 Convert metagenomSeq object with normolized count into dataframe  
css.otu.all <- data.frame(otu_table(metagenomeSeq::MRcounts(mg.ps1.1, norm = TRUE, log = TRUE), taxa_are_rows = TRUE))

#6.1.6 Subset only significant taxa 
css.otu <- t(css.otu.all[colnames(rfdata) %in% rownames(pval.rf.perm.s1), ])

#6.1.7 Add taxa collnames 
colnames(css.otu) <- rownames(pval.rf.perm.s1)

#6.1.7 Add metadata 
css.otu.comb <- cbind(css.otu, ps1.1@sam_data[ ,c("Month", "Year", "CowID", "Status", "Shading")])

#6.1.7 Add adjusted collumn 
css.otu.comb$Year_month <- paste(css.otu.comb$Year, css.otu.comb$Month, sep = "_")

#6.1.8 Melt into long format 
css.otu.comb.m <- reshape2::melt(css.otu.comb)

#6.1.9 Correct spelling of the Shedding column name 
colnames(css.otu.comb.m)[5] <- "Shedding"

rm(list = c("mg.ps1.1", "p", "css.otu.comb"))

#6.2 Plot normalized abundace
##############################

#6.2.1 Make an empty list to store plots 
plots.all <- list()

#6.2.2 Make run loop to create box and time serious plots 
for (taxa in levels(css.otu.comb.m$variable)) {
  
  # Subset taxa of interest 
  data.gg <- css.otu.comb.m[css.otu.comb.m$variable %in% taxa, ]
  
  # time seriouse plot 
  time.gp.tax <- ggplot(data.gg, aes(y=value, x=Year_month, color = CowID, 
                                     shape = Shedding, group = CowID, linetype=Status)) + 
                 theme_bw() + 
                 geom_point(size=2, alpha = 0.7, stroke = 1) + 
                 scale_color_brewer(palette = "Dark2") + 
                 geom_line(size=1) + 
                 theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
                 ggtitle(taxa) + 
                 scale_shape_manual(values=c(15,16,17,18))
  
  # box plot 
  box.gp.tax <- ggplot(data.gg,  
                aes(y=value, x=Status)) + 
                theme_bw() + 
                geom_boxplot() + 
                geom_jitter(data = css.otu.comb.m[css.otu.comb.m$variable %in% taxa, ], 
                aes(y=value, x=Status, color = CowID, shape = Shedding), 
                width = 0.2, size = 2, stroke = 1,alpha = 0.7) +
                scale_color_brewer(palette = "Dark2") + 
                ggtitle(taxa) + 
                scale_shape_manual(values=c(15,16,17,18))
  # combine plots 
  time.box.comb <- cowplot::plot_grid(time.gp.tax, box.gp.tax, ncol = 2, rel_widths = c(4/6, 2/6))
  
  # add combine plots to list 
  plots.all[[taxa]] <- time.box.comb
}

# 6.2.3 Save plots 
pdf("output/plots/Supplementary_plots1.pdf", width = 12, height = 5)
invisible(lapply(plots.all, print))
dev.off()


#6.3 Plot relative abundance in the same fromat 
###############################################

#6.3.1 Subset singnificantly contributing taxa from relative abundace table (rfdata)
relat.otu <- rfdata[, colnames(rfdata) %in% rownames(pval.rf.perm.s1) ]

#6.3.2 Add value relative abundance to dataframe created at step 6.2 
css.otu.comb.m$relat <- reshape::melt(relat.otu)$value

#6.3.3 Make an empty list for plot objects 
plots.all <- list()

#6.3.4 Plot objects in a loop 
for (taxa in levels(css.otu.comb.m$variable)) {
  
  # Subset a single taxa 
  data.gg <- css.otu.comb.m[css.otu.comb.m$variable %in% taxa, ]
  
  # Plot time serious 
  time.gp.tax <- ggplot(data.gg, aes(y=relat, x=Year_month, color = CowID, 
                                     shape = Shedding, group = CowID, linetype=Status)) + 
                  theme_bw() + 
                  geom_point(size=2, alpha = 0.7, stroke = 1) + 
                  scale_color_brewer(palette = "Dark2") + 
                  geom_line(size=1) + 
                  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
                  ggtitle(taxa) + 
                  scale_shape_manual(values=c(15,16,17,18))
  
  # Plot boxplot 
  box.gp.tax <- ggplot(data.gg,  
                       aes(y=relat, x=Status)) + 
                      theme_bw() + 
                      geom_boxplot() + 
                      geom_jitter(data = css.otu.comb.m[css.otu.comb.m$variable %in% taxa, ], 
                      aes(y=value, x=Status, color = CowID, shape = Shedding), 
                      width = 0.2, size = 2, stroke = 1,alpha = 0.7) +
                      scale_color_brewer(palette = "Dark2") + 
                      ggtitle(taxa) + 
                      scale_shape_manual(values=c(15,16,17,18))
  
  # Combine plots 
  time.box.comb <- cowplot::plot_grid(time.gp.tax, box.gp.tax, ncol = 2, rel_widths = c(4/6, 2/6))
  
  # Add combine plot in the list 
  plots.all[[taxa]] <- time.box.comb
}

#6.3.5 Save objects 
pdf("output/plots/Supplementary_plots1_relative.pdf", width = 12, height = 5)
invisible(lapply(plots.all, print))
dev.off()



#7. LEFSe analisis (external, run on galaxy server)
###################################################

#7.1 Prepare data for LEFSe analysis 
lefse.data <- data.frame(t(rfdata[, !colnames(rfdata) %in% "Status"]))

lefse.data <- rbind(as.character(rfdata[, colnames(rfdata) %in% "Status"]), lefse.data)

#7.1.1 Write data. Data should be saved via ecxel to have the correct format for LEFSe
write.table(lefse.data, "output/processing/LEFse_input.txt", col.names=FALSE, sep = "\t")


#7.2 Read output data from LEFSe (defult parameteres used for biomarkers identification)
lefse.out <- read.delim("output/LEFSe_external/Galaxy9-[B)_LDA_Effect_Size_(LEfSe)_on_data_8].lefse_internal_res", 
              header = FALSE)

lefse.out$ASV.ID <- gsub("^.*_", "", lefse.out$V1)


#7.3 Add extra data 
#7.3.1 Create a named (taxa ID as names) vector that contain sequences 
all.seq <- rownames(data.frame(ps1.1@tax_table))

names(all.seq) <- colnames(rfdata)[1:(ncol(rfdata)-1)]

#7.3.2 Convert the object into a dataframe that will be combined with lefse output 
all.seq.df <- data.frame(cbind(all.seq, names(all.seq), gsub("^.*[.]", "", names(all.seq))))

colnames(all.seq.df) <- c("Seq", "ID", "ASV.ID")

#7.3.3 Combine Dataframes 
lefse.out.1 <- left_join(lefse.out, all.seq.df, by="ASV.ID")



#8.0 Bult general phylogenetic tree 
####################################

#8.1 Use named vector with sequenses 

#8.2 Run alligment 
#8.2.1 Load library 
library(msa)

#8.2.2 Aligment 
seq.aligment <- msa(all.seq, method="ClustalW", type="dna", order="input")


#8.2.3 Save alligment file 
#save(seq.aligment, file = "output/processing/seq_aligment.Rdata")
detach("package:msa", unload=TRUE)


#8.3 Bult phylogenetic tree (not rooted)
#8.3.1 Convert aligment object to phyDat format 
phang.align <- as.phyDat(seq.aligment, type="DNA", names=getSequence(seq.aligment))

#8.3.2 Culculate distances between sequenses 
dm <- dist.ml(phang.align)

#8.3.3 Bult first nebour joined tree 
treeNJ <- phangorn::NJ(dm) # Note, tip order != sequence order

#8.3.4 Fit the tree 
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fit, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

#8.3.5 Save output object 
#save(fitGTR, file = "output/processing/fitGTR.Rdata")


#8.4 Visualize phylogenetic tree 

#8.4.1 Extract taxonomy from a corresponding phyloseq object 
tax.tree <- data.frame(ps1.1@tax_table)

#8.4.2 Add names of the tip lables from the tree 
tax.tree$tipID <- fitGTR$tree$tip.label

#8.4.3 Make a list with with members of certain phylums arranged in vectors
color.tree <- list()
for (p in levels(tax.tree$Phylum)) { 
  color.tree[[p]]<- tax.tree$tipID[tax.tree$Phylum %in% p]
  }

#8.4.4 Arrange list by number of enteries in each vector
color.tree <- color.tree[order(sapply(color.tree,length),decreasing=T)]

#8.4.5 Make color pallet for phylum 
palette.phylum  <- RColorBrewer::brewer.pal(8,"Dark2") 
palette.phylum.f  <- c( palette.phylum, rep("gray80", 50)) 

#8.4.6 Add names phylum names as color names 
names(palette.phylum.f) <- c(names(color.tree), "0")

#8.4.7 Root the tree at midpoint to plot it incirular fation 
root.t <- phangorn::midpoint(fitGTR$tree)

#8.4.8 Make initial tree plot 
p.tree <- ggtree::ggtree(root.t, layout="circular", branch.length='none')

#8.4.9 Group ASV by phylum 
p.gr.phylum <- groupOTU(p.tree, color.tree, 'Phylum') 

#8.4.10 Rearange levels within color groups (phylum) 
p.gr.phylum$data$Phylum <- factor(p.gr.phylum$data$Phylum, 
                                  levels = c(names(color.tree), "0"))

#8.4.11 Plot the plot 
p.gr.phylum + aes(color=Phylum) +
  theme(legend.position="right")  + scale_colour_manual(values=palette.phylum.f) 

ggsave("output/plots/phy_tree_all.pdf")
ggsave("output/plots/phy_tree_all.png", dpi = 300)





#8.5 Select and plot tree only with biomarkers 
###############################################

source("src/sup_f_tree_plot.R")
#8.5.1 Subset the tree 
biom.sig.tree <- drop.tip(fitGTR$tree, fitGTR$tree$tip.label[!fitGTR$tree$tip.label %in% rownames(pval.rf.perm.s1)])

#8.5.1 Root trimmed tree at midpoint 
biom.sig.tree.r <- midpoint(biom.sig.tree)

#8.5.2 Prepare coloring data 
# a. Extract taxomic information for sig taxa
tax.tree.sig <- tax.tree[tax.tree$tipID %in% biom.sig.tree.r$tip.label, ]

# b. Make a list with with members of certain phylums arranged in vectors
color.tree.sig <- list()
for (p in levels(tax.tree.sig$Phylum)) { 
  if (length(tax.tree.sig$tipID[tax.tree.sig$Phylum %in% p]) > 0) {
    color.tree.sig[[p]]<- tax.tree.sig$tipID[tax.tree.sig$Phylum %in% p]
  }
}

#8.5.3 Make initial tree plot 
p.tree.sig <- ggtree(biom.sig.tree.r, size=1)

#8.5.4 Group ASV by phylum 
p.gr.phylum.sig <- groupOTU(p.tree.sig, color.tree.sig, 'Phylum') 

#8.5.5 Rearange levels within color groups (phylum) 
p.gr.phylum.sig$data$Phylum <- factor(p.gr.phylum.sig$data$Phylum, 
                        levels = c(names(color.tree.sig[order(sapply(color.tree.sig,length),decreasing=T)]), "0"))

#8.5.6 Plot the plot 
tree.p <- p.gr.phylum.sig + aes(color=Phylum) +
  theme(legend.position="right")  + scale_colour_manual(values=palette.phylum.f) + 
  geom_tiplab(size=4, align=TRUE, linesize=.5) + xlim(0, 1.2) + scale_y_tree() + 
  theme(legend.position = "bottom")


#8.6 Prepare data to add additional information 
#8.6.1 Add data abbout MDA 
rf.ant <- data.frame(cbind(rownames(pval.rf.perm.s1), pval.rf.perm.s1$MeanDecreaseAccuracy.scaled))

rf.ant$MDA.hl <- pval.rf.perm.s1$High.scaled -  pval.rf.perm.s1$Low.scaled

colnames(rf.ant) <- c("label", "MDA.score", "High.Low")

#8.6.2 Add data abbout Lefse score 
lf.ant <- lefse.out.1[lefse.out.1$ID %in% rownames(pval.rf.perm.s1), c("V4", "V3", "ID")]

colnames(lf.ant) <- c("LefSe.score", "Indicate", "label")


#8.6.3 Add data about avarage abundace in low shedders 
abund.sig.l <- rfdata[ps1.1@sam_data$Status %in% "Low", colnames(rfdata) %in% rownames(pval.rf.perm.s1)]

mean.ab.sig.l <- cbind(data.frame(rowMeans(as.matrix(t(abund.sig.l)))*100), colnames(abund.sig.l))

colnames(mean.ab.sig.l) <- c("Mean.abund", "label")

mean.ab.sig.l$Group <- "Low"



#8.6.4 Add data about avarage abundace in high shedders 
abund.sig.h <- rfdata[ps1.1@sam_data$Status %in% "High", colnames(rfdata) %in% rownames(pval.rf.perm.s1)]

mean.ab.sig.h <- cbind(data.frame(rowMeans(as.matrix(t(abund.sig.h)))*100), colnames(abund.sig.h))

colnames(mean.ab.sig.h) <- c("Mean.abund", "label")

mean.ab.sig.h$Group <- "High"



#8.6.5 Add data about prevalence in low shedders 
prev.sig.l <- abund.sig.l

prev.sig.l[prev.sig.l > 0] <- 1

prev.t.sig.l <- cbind(data.frame(round(colSums(prev.sig.l)/nrow(prev.sig.l)*100, 2)), colnames(prev.sig.l))

colnames(prev.t.sig.l) <- c("Prevalence", "label") 

prev.t.sig.l$Group <- "Low"


#8.6.6 Add data about prevalence in high shedders
prev.sig.h <- abund.sig.h

prev.sig.h[prev.sig.h > 0] <- 1

prev.t.sig.h <- cbind(data.frame(round(colSums(prev.sig.h)/nrow(prev.sig.h)*100, 2)), colnames(prev.sig.h))

colnames(prev.t.sig.h) <- c("Prevalence", "label") 

prev.t.sig.h$Group <- "High"



#8.7 Arrange and plot addtional iformation 
##########################################

  # Extract tip order 
tip.order <- filter(tree.p, isTip) %>% select(c(label, y))


#8.7.1 Ad information about MDA score 
library("viridis")    

rf.ant.o <- left_join(rf.ant, tip.order, by='label')

rf.ant.p <- ggtreeplot(tree.p, rf.ant.o , aes(y=MDA.score, width=0.8), flip=TRUE) +
  geom_col(aes(fill=High.Low)) +
  coord_flip() + 
  theme_tree2() + 
  geom_vline(xintercept = seq(1.5, 70, 1), color = "grey92") + 
  clean_x() +
  scale_fill_viridis(option = "D") +
  theme(legend.position = "bottom")


#8.7.2 Add information about lefse score 

lf.ant.o <- left_join(lf.ant, tip.order, by='label')

lf.ant.p <- ggtreeplot(tree.p, lf.ant.o , aes(y=LefSe.score, width=0.8), flip=TRUE) +
  geom_col(aes(fill=Indicate)) +
  coord_flip() + 
  theme_tree2() + 
  geom_vline(xintercept = seq(1.5, 70, 1), color = "grey92") + 
  clean_x() +
  #scale_fill_viridis(option = "D") +
  theme(legend.position = "bottom")


#8.7.3 Add information about mean abundunce in low group 
mean.ab.sig.l.o <- left_join(mean.ab.sig.l, tip.order, by='label')

mean.ab.sig.l.p <- ggtreeplot(tree.p, mean.ab.sig.l.o , aes(y=Mean.abund, width=0.8), flip=TRUE) +
  geom_col(aes(fill=Group)) +
  coord_flip() + 
  theme_tree2() + 
  geom_vline(xintercept = seq(1.5, 70, 1), color = "grey92") + 
  clean_x() +
  theme(legend.position = "bottom") + scale_fill_manual(values = "Steelblue")


#8.7.4 Add information about mean prevalece in low group 
prev.t.sig.l.o <- left_join(prev.t.sig.l, tip.order, by='label')

prev.t.sig.l.p <- ggtreeplot(tree.p, prev.t.sig.l.o, aes(y=Prevalence, width=0.8), flip=TRUE) +
  geom_col(aes(fill=Group)) +
  coord_flip() + 
  theme_tree2() + 
  geom_vline(xintercept = seq(1.5, 70, 1), color = "grey92") + 
  clean_x() +
  theme(legend.position = "bottom") + scale_fill_manual(values = "Steelblue")


#8.7.5 Add information about mean abundunce in high group
mean.ab.sig.h.o <- left_join(mean.ab.sig.h, tip.order, by='label')

mean.ab.sig.h.p <- ggtreeplot(tree.p, mean.ab.sig.h.o , aes(y=Mean.abund, width = 0.8), flip=TRUE) +
  geom_col(aes(fill=Group)) +
  coord_flip() + 
  theme_tree2() + 
  geom_vline(xintercept = seq(1.5, 70, 1), color = "grey92") + 
  clean_x() +
  theme(legend.position = "bottom") + 
  scale_fill_manual(values = "Red") + expand_limits(y = 0.6)


#8.7.6 Add information about mean prevalence in high group
prev.t.sig.h.o <- left_join(prev.t.sig.h, tip.order, by='label')

prev.t.sig.h.p <- ggtreeplot(tree.p, prev.t.sig.h.o, aes(y=Prevalence, width=0.8), flip=TRUE) +
  geom_col(aes(fill=Group)) +
  coord_flip() + 
  theme_tree2() + 
  geom_vline(xintercept = seq(1.5, 70, 1), color = "grey92") + 
  clean_x() +
  theme(legend.position = "bottom") + scale_fill_manual(values = "Red")


#8.7.7 Plot plot data togather 
cowplot::plot_grid(tree.p, 
                   rf.ant.p, 
                   lf.ant.p, 
                   mean.ab.sig.l.p,
                   prev.t.sig.l.p,  
                   mean.ab.sig.h.p , 
                   prev.t.sig.h.p, 
                      ncol=7, align='h', 
                      labels=LETTERS[1:7], 
                      rel_widths = c(1, .25, .25, .15, .15, .15, .15))

# Save resulted plot 
ggsave("output/plots/biomarkers.pdf", height = 12, width = 16)

