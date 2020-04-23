#
#######################################################
# Random forest classification and results visualiation 
# on Genus level 
#######################################################


# 0 Prepare enviorment 
########################

# 0.1 Set working directory 

setwd("C://Rprojects/Microbiota_MAP/")
dir.create("output/rf_objects/genus_level")

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



################################################
# 1 Preapareation of input data for RF analysis 
################################################

# 1.1 Import filtered phyloseq object 
ps1.1 <- readRDS("output/processing/ps1.1.rds")

# 1.2 Summarize taxa to genus level
ps1.1.genus <- tax_glom(ps1.1, "Genus", NArm = FALSE)
  
# 1.2 Convertion of the data to proportional data 
ps1.1.genus.prop <- transform_sample_counts(ps1.1.genus, function(otu) otu/sum(otu))

# 1.3 Extraction of OTU tabel 
rfdata.genus <- data.frame(ps1.1.genus.prop@otu_table@.Data)

# 1.4 Format genus names.   
#     Names of unssigned genera(UG) will composed from the last assinged taxonomic level and UG and count number 
# 1.4.1 Extract taxonomic table 
ps1.1tax.genus <- data.frame(ps1.1.genus.prop@tax_table)

# 1.4.2 Make a vector containing last ide 
ug.vect <- c()
for (r in 1:nrow(ps1.1tax.genus)) {
  ind.r <- ps1.1tax.genus[r, ]
  ug.vect <- c(ug.vect, tail(ind.r[!is.na(ind.r)], 1))
}
rm(list = c("r", "ind.r"))

# 1.5 Add last assigned taxonomic level and ASV number as column names for OTU table 

colnames(rfdata.genus) <- paste(gsub("-|/", "", ug.vect), 1:nrow(ps1.1tax.genus), sep = ".")
rm(list = c("ug.vect", "ps1.1tax"))  

# 1.6 Add to rf data column with status 
rfdata.genus$Status <- ps1.1.genus.prop@sam_data$Status



################################################
# 2. Optimization of parameters for RF model 
################################################

#2.1 Optimiztion of nubber of features used per split used for with different number of trees. 
##################################################################################

#2.1.1 list number of used trees 
ntrees <- c(1000, 5000, 10000, 15000, 20000)

#2.1.2 Create an empty list 
mtry.out2.g <- list()

#2.1.3 Run the adjustment of mtry parrametar with various number of trees 
#Set seed 
set.seed(4574)
sqrt(ncol(rfdata.genus))

for (tr in as.character(ntrees)) {
  mtry.out1.g <- foreach(i=1:8) %dopar% {  randomForest::tuneRF(rfdata.genus[, -which(names(rfdata.genus) %in% "Status")],
                                                                rfdata.genus[, which(names(rfdata.genus) %in% "Status")],
                                                              13, ntreeTry=as.numeric(tr), stepFactor=2, improve=0.05,
                                                              trace=TRUE, plot=TRUE, doBest=FALSE)
  }
  mtry.out2.g[[tr]] <-  mtry.out1.g
}
rm(list = c("tr", "i", "mtry.out1.g"))

#2.1.4 Save resulted object 
save(mtry.out2.g, file="output/rf_objects/genus_level/mtry_out2_g.Rdata")

#2.1.5 Convert the object in to dataframe 
mtry.out2.g.df <- c()
for (tr in as.character(ntrees)) {
  for (i in 1:8) {
    one.try <- data.frame(mtry.out2.g[[tr]][[i]])
    one.try$ntrees <- tr
    one.try$run <- i
    
    mtry.out2.g.df <- rbind(mtry.out2.g.df, one.try)
  }
}
rm(list = c("tr", "i", "one.try"))
gc()

#2.1.6 Adjust dataframe for correct plotting 
mtry.out2.g.df$OOBError.round <- round(as.numeric(mtry.out2.g.df$OOBError), 2)

mtry.out2.g.df$mtry <- factor(mtry.out2.g.df$mtry, 
                            levels = as.character(unique(mtry.out2.g.df$mtry)[order(unique(mtry.out2.g.df$mtry))]))
mtry.out2.g.df$ntrees <- factor(mtry.out2.g.df$ntrees, 
                              levels = c("1000", "5000", "10000", "15000", "20000"))


#2.1.7 Plot results 
plot.mtr.ntr.g <- ggplot(data = mtry.out2.g.df, aes(x=mtry, y=OOBError.round, group=as.factor(run), color=as.factor(run))) + 
  geom_jitter(size=3, alpha = 0.5, width = 0.05) + 
  geom_line(alpha = 0.7) +
  facet_grid(ntrees ~., scales = "free") + 
  theme_bw()

#2.1.8 Save the plot 
ggsave("output/processing/ntree_mtry_asv_genus.pdf", plot = plot.mtr.ntr.g)
rm(list = c("plot.mtr.ntr", "mtry.out2.df"))


#2.2 Buld RF models in 8 folds using 7 per split (best parameters)
#2.2.1 Buld 8 RF models with each number of trees and mtry=176
ntrees.out.g <-list()
for(nt in as.character(ntrees)) {
  rf.rep <- foreach(i=1:8) %dopar% {  randomForest::randomForest(Status ~ ., 
                                                                 data = rfdata.genus, 
                                                                 importance=TRUE,
                                                                 proximity=TRUE, 
                                                                 mtry=7,
                                                                 ntree=as.numeric(nt))
  }
  ntrees.out.g[[nt]] <- rf.rep
  gc()
}


#2.2.2 Remove temp loop objects 
rm(list=c("nt", "ntrees", "rf.rep"))

#2.2.3 Save resulted object 
save(ntrees.out.g,
     file = "output/rf_objects/genus_level/ntrees_out_g.RData")

#2.2.5 Extract information about class errors 
rf.ntree <- c()
for (rf1 in names(ntrees.out.g)) {
  for (rf2 in 1:8) {
    rf.out <- data.frame(ntrees.out.g[[rf1]][[rf2]]$confusion)
    rf.out$ntree <- rf1
    rf.out$nmtry <- "splist@7"
    rf.out$Status <- rownames(rf.out)
    
    rf.ntree <- rbind(rf.ntree, rf.out)
  }
}
rm(list=c("rf1", "rf1.out", "rf2"))
gc()

#2.2.6 Convert informaiton about class errors into dataframe 
rf.ntree.df <- data.frame(rf.ntree)
rf.ntree.df$ntree <- factor(rf.ntree.df$ntree, 
                            levels = c("1000", "5000", "10000", "15000", "20000"))
#2.2.7 Plot class errors 
best.tree.split <- ggplot(rf.ntree.df, aes(y= as.numeric(as.character(class.error)),
                                           x=ntree, color = Status)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1, height = 0.001, size = 3, alpha=0.75) + 
  theme_bw() + 
  facet_grid(~nmtry) +
  ylab("Error rate")

#2.2.8 Save the plot 
ggsave("output/processing/best_tree_split_genus.pdf", plot = best.tree.split)                       
rm(list=c("rf1", "rf1.out", "rf2", "rf.out2")) 


#3. Influence of the dataset size on accuracy of classification 
###############################################################

#3.1 Define function that will subest defined presetage of samples from High and Low 
#    shading groups 
ssub.function <- function(ssub) {
  rfdata.genus.high <- rfdata.genus[rfdata.genus$Status %in% "High", ]
  rfdata.genus.high.ssub <- rfdata.genus.high[sample(nrow(rfdata.genus.high), 
                                         ssub*nrow(rfdata.genus.high), 
                                         replace = FALSE), ]
  rfdata.genus.low <- rfdata.genus[rfdata.genus$Status %in% "Low", ]
  rfdata.genus.low.ssub <- rfdata.genus.low[sample(nrow(rfdata.genus.low), 
                                       ssub*nrow(rfdata.genus.low), 
                                       replace = FALSE), ]
  rfdata.genus.ssub <- rbind(rfdata.genus.high.ssub , rfdata.genus.low.ssub)
  return(rfdata.genus.ssub)
}

#3.2 Run 16 alteration of the RF with different number of subsetted samples
#3.2.1 Define eviorment and variables 
ssub.out <-list()
ssub.v <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4)
set.seed(24735)

#3.2.2 Loop itself 
for(ssub in ssub.v) {
  rf.rep <- foreach(i=1:16) %dopar% {randomForest::randomForest(Status ~ ., 
                                                                data = ssub.function(ssub), 
                                                                importance=TRUE,
                                                                proximity=TRUE, 
                                                                mtry=7,
                                                                ntree=10000)
  }
  ssub.out[[as.character(ssub)]] <- rf.rep
  gc()
}

#3.2.3 Save resulted object 
save(ssub.out, file = "output/rf_objects/ssub_out_genus.RData")

#3.3 Extract data from the resulted object 
ssub.data <- c()
for (rf1 in names(ssub.out)) {
  for (rf2 in 1:16) {
    rf.out <- data.frame(ssub.out[[rf1]][[rf2]]$confusion)
    rf.out$ntree <- rf1
    rf.out$Status <- rownames(rf.out)
    
    ssub.data <- rbind(ssub.data, rf.out)
  }
}
rm(list=c("rf1", "rf1.out", "rf2", "rf.out2"))


#3.4 Plot class errors per group with fitted line (Locally Weighted Regression)
ssub.plot <- ggplot(ssub.data, aes(y= as.numeric(as.character(class.error)),
                                   x=as.numeric(ntree), color = Status)) + 
  geom_jitter(width = 0.01, height = 0.001, size = 2, alpha=0.75) + 
  stat_smooth(method = "loess", formula = y ~ x, size = 1) + 
  theme_bw() + 
  ylab("Error rate") + 
  xlab("Proportion of the data")

#3.5 Save the plot 
ggsave("output/processing/data_size_ifluence_genus.pdf", plot = ssub.plot) 


