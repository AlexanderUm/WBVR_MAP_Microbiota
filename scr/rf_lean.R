# Lean RF model
###############

# This script will ascees: 
# 1. Accuracy of the RF model that used only ASV classified as stasitsitcall significant descriminators (biomarkers) 
# 2. Accuracy of the model when genera that coplite genera containing (biomarker)

#I Accurary of the RF model 

#0. Set up enviorment 
# 0.1 Set working directory 

setwd("C://Rprojects/Microbiota_MAP/")

# 0.2 Load libraries 
library(randomForest)
library(ggplot2)


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


source("src/rf_funcitons.R")
#1. Import and subset the rf object 
load("output/processing/rfdata.Rdata")

load("output/processing/pval_rf_perm.Rdata")

rfdata.biom <- rfdata[ , colnames(rfdata) %in% c(rownames(pval.rf.perm.s1), "Status")]



#I. Constract rf using only significant taxa (RF identified )
##################################################

#I. 1 number of trees and mtry parameter optimization 
ntrees <- c(1000, 5000, 10000, 15000, 20000)

nt.mtr.rfbiom <- TreeMtryPlot(data=rfdata.biom, ntrees = ntrees, class_colum = "Status", ntimes = 7)

ggsave("output/plots/ntree_mtr_rfbiom.pdf", plot = nt.mtr.rfbiom, width = 7, height = 7)
ggsave("output/plots/ntree_mtr_rfbiom.png", plot = nt.mtr.rfbiom, dpi = 300, width = 7, height = 7)


#I. 2 Class error with best splits 
accur.bs.rfbiom <- BestSplitAccuracy(data = rfdata.biom, ntrees = ntrees, mtry = c(8, 16), class_colum = "Status", ntimes = 7)

ggsave("output/plots/accuracy_bs.pdf", plot = accur.bs.rfbiom, width = 7, height = 5.5)
ggsave("output/plots/accuracy_bs.png", plot = accur.bs.rfbiom, dpi = 300, width = 7, height = 5.5)


#I. 3 Influence of datasize on accuracy of classification 
sub.percent <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4)

data.s.infl.rfbiom <- PlotSizeInfluece(data = rfdata.biom, 
                                        sub.percent = sub.percent, 
                                        class_colum = "Status",  
                                        mtr = 8, ntree = 15000, nrep = 16)

ggsave("output/plots/data_size_infl_rfbiom.pdf", plot = data.s.infl.rfbiom, width = 7, height = 4)
ggsave("output/plots/data_size_infl_rfbiom.png", plot = data.s.infl.rfbiom, dpi = 300,  width = 7, height = 4)



#II. Using taxa identified by LefSe as biomarkers to build RF modles 
#######################################################################
#II. 0 Import data from LefSe and perepare data for analysis 
lefse.out <- read.delim("output/LEFSe_external/Galaxy9-[B)_LDA_Effect_Size_(LEfSe)_on_data_8].lefse_internal_res", 
                        header = FALSE) 

rfdata2 <- rfdata

colnames(rfdata2) <- gsub("[.]", "_", colnames(rfdata2))

rfdata2 <- rfdata2[colnames(rfdata2) %in% lefse.out$V1[!is.na(lefse.out$V4)]]

rfdata2$Status <- rfdata$Status


#II. 1 number of trees and mtry parameter optimization 
nt.mtr.lefse <- TreeMtryPlot(data=rfdata2, ntrees = ntrees, class_colum = "Status", ntimes = 7)

ggsave("output/plots/ntree_mtr_lefse.pdf", plot = nt.mtr.lefse, width = 7, height = 7)
ggsave("output/plots/ntree_mtr_lefse.png", plot = nt.mtr.lefse, width = 7, height = 7, dpi = 300)


#II. 2 Influence of datasize on accuracy of classification 
accur.bs.lefse <- BestSplitAccuracy(data = rfdata2, ntrees = ntrees, mtry = c(15, 30), class_colum = "Status", ntimes = 7)

ggsave(filename = "output/plots/accuracy_bs_lefse.pdf", plot = accur.bs.lefse, width = 7, height = 5.5)
ggsave(filename = "output/plots/accuracy_bs_lefse.png", plot = accur.bs.lefse, width = 7, height = 5.5, dpi=300)


#II. 3 Influence of datasize on accuracy of classification
data.s.infl.lefse <- PlotSizeInfluece(data = rfdata2, 
                                        sub.percent = sub.percent, 
                                        class_colum = "Status",  
                                        mtr = 30, ntree = 15000, nrep = 16)

ggsave("output/plots/data_size_infl_lefse.pdf", plot = data.s.infl.lefse, width = 7, height = 4)
ggsave("output/plots/data_size_infl_lefse.png", plot = data.s.infl.lefse, dpi = 300, width = 7, height = 4)


#III. Use as predictior genera where significan contributing ASV were detected 
##############################################################################

# 1.1 Import filtered phyloseq object 
ps1.1 <- readRDS("output/processing/ps1.1.rds")


#2.1 Glom timmed ps at genus level 
ps1.1.genus <- phyloseq::tax_glom(ps1.1, "Genus", NArm = FALSE)

ps1.1.g.prop <- phyloseq::transform_sample_counts(ps1.1.genus, function(otu) otu/sum(otu))


ps.biom <- prune_samples(colnames(rfdata)["Status"] %in% rownames(pval.rf.perm.s1), ps1.1) 


tax.t.genus <- data.frame(ps1.1.genus@tax_table)


for (lvl in colnames(tax.t.genus)) {
  
  tax.t.genus <- data.frame(ps1.1.p.b@tax_table)
  
  ps1.1.p.b <-  prune_taxa(tax.t.genus[, lvl] %in% tax.t.genus.biom[, lvl], ps1.1.p.b)
  
}




#
# 1.2 Summarize taxa to genus level
ps1.1.genus <- phyloseq::tax_glom(ps1.1, "Genus", NArm = FALSE)

# 1.3 Extract both tax tables 


ps1.1.p.b2 <- prune_taxa(as.vector(!is.na(ps1.1.p.b@tax_table[,"Genus"])), ps1.1.p.b) 

rf.ps1.1pb <- data.frame(ps1.1.p.b2@otu_table)
colnames(rf.ps1.1pb) <- as.vector(ps1.1.p.b2@tax_table[, "Genus"])
rf.ps1.1pb$Status <- ps1.1.p.b@sam_data$Status
  


rf.ps1.1pb <- data.frame(ps1.1.p.b@otu_table)
colnames(rf.ps1.1pb) <- paste(as.vector(ps1.1.p.b@tax_table[, "Genus"]), 1:ncol(rf.ps1.1pb), sep = "_G")
rf.ps1.1pb$Status <- ps1.1.p.b@sam_data$Statu

ntrees <- c(1000, 5000, 10000, 15000, 20000)

TreeMtryPlot(data=rf.ps1.1pb, ntrees = ntrees, class_colum = "Status", ntimes = 7)



BestSplitAccuracy(data = rf.ps1.1pb, ntrees = ntrees, mtry = c(4, 7), class_colum = "Status", ntimes = 7)


sub.percent <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4)

PlotSizeInfluece(data = rf.ps1.1pb, 
                 sub.percent = sub.percent, 
                 class_colum = "Status",  
                 mtr = 4, ntree = 5000, nrep = 16)







