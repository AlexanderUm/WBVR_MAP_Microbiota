#Classification of new samples 

#Set enviorment 
################

set.seed(2947)
setwd("C://Rprojects/Microbiota_MAP_II_20200417/")

library(randomForest)
library(tidyverse)
library(phyloseq)

#1. Retrive objects 
################



#RF models generted with the first batch 
load("resource/RF_models.R")



#1.1 Investigate data

#First phyloseq PCoA
################

load("resource/phyloseq0_firstRun.RData")
ps0.I <- ps
ps0.I.prop <- transform_sample_counts(ps0.I, function(otu) otu/sum(otu))


colPCoA <- RColorBrewer::brewer.pal(4, "Dark2")
colPCoA <- c(colPCoA, "red", "red", replicate(n=10, "grey"))

GP.ord <- ordinate(ps0.I.prop, "MDS", "bray")
p1 = plot_ordination(ps0.I.prop, GP.ord, type="samples", color = "Shading", 
                     title="Samples")

PCoA.I <- p1 + geom_point(size=3) + theme_bw() + scale_color_manual(values = colPCoA)

ggsave("output/test_plots/PCoA_All_I.pdf", PCoA.I)
ggsave("output/test_plots/PCoA_All_I.png", PCoA.I, dpi = 300)


#Second phyloseq 
#Phyloseq of the second run
load("output/dada2/phyloseq0.RData")

ps0.II <- ps 

ps0.II.prop <- transform_sample_counts(ps0.II, function(otu) otu/sum(otu))

psII.meta <- ps0.II.prop@sam_data
psII.meta$value <- as.character(psII.meta$value)
psII.meta$value[is.na(psII.meta$value)] <- psII.meta$Sample.name[is.na(psII.meta$value)]
psII.meta$value[psII.meta$value %in% "neg"] <- "Blank"
ps0.II.prop@sam_data <- psII.meta

colPCoA2 <- RColorBrewer::brewer.pal(4, "Dark2")
colPCoA2 <- c(colPCoA2, "red", replicate(n=10, "grey"))

GP.ord <- ordinate(ps0.II.prop, "MDS", "bray")
p1 = plot_ordination(ps0.II.prop, GP.ord, type="samples", color = "value", 
                     title="Samples")

PCoA.II <- p1 + geom_point(size=3) + theme_bw() + scale_color_manual(values = colPCoA2)

ggsave("output/test_plots/PCoA_All_II.pdf", PCoA.II)
ggsave("output/test_plots/PCoA_All_II.png", PCoA.II, dpi = 300)


#Classify samples using developed moddels 
##########################################

#Extract data from phyloseq and rf.lean.model
rfdataII <- data.frame(ps0.II.prop@otu_table@.Data)

rf.lean.taxa <- rownames(rf.lean.model$importance)

#Make objects for the loop 
overl.lean.tax <- list()

overl.full.tax <- list()

samp.sammary <- data.frame(ID = character(), 
                   TotalTax = numeric(), 
                   Overl_lean = numeric(), 
                   Overl_full = numeric(), 
                   Class_lean = character(), 
                   Class_full = character(), 
                   Class_indRF = character(),
                   OOBE_indRF = character(), 
                   Low_indRF = character(), 
                   High_indRF = character(),
                   Mtry_indRF = character(), 
                   Sheddiing_score = numeric(), 
                   stringsAsFactors = FALSE)
mtry.steps <- list()

rf.models <-list()

#Loop: subset indididual samples; find a number of nonzero samples per taxa
#      find overlap with taxa in lean and full models; build RF using overlaping taxa and 
#      dataset I; fill out the sammary table 

for (i in 1:nrow(rfdataII)) {
  #remove taxa that are not present 
  samp0 <- rfdataII[i, ]
   
  samp <- samp[, !samp == 0]
  
  samp.sammary[i, "TotalTax"] <- ncol(samp)
  
  #find overlap with taxa in lean and full rf models 
  overl.lean.tax[[i]] <- intersect(rf.lean.taxa, colnames(samp))
  
  samp.sammary[i, "Overl_lean"] <- length(intersect(rf.lean.taxa, colnames(samp)))
  
  overl.full.tax[[i]] <- intersect(rf.full.taxa, colnames(samp))
  
  samp.sammary[i, "Overl_full"] <- length(intersect(rf.full.taxa, colnames(samp)))
  
  #Build the RF model based on samples intersect
  rfdata.samp <- rfdata.seqnames[, c(intersect(rf.full.taxa, colnames(samp)), "Status") ]
 
      #adjust number of trees used in a basket 
    vec.tr <- c(round(sqrt(ncol(samp)), 0), 
              round(ncol(samp)*0.1, 0), 
              round(ncol(samp)*0.03, 0), 
              round(ncol(samp)*0.15, 0)
              ) 
    
     bestmtry.all <- c()
     
    for (tr in vec.tr) {
      bestmtry <- tuneRF(rfdata.samp[, !colnames(rfdata.samp) %in% "Status"], 
                       rfdata.samp[,"Status"],
                       mtryStart = tr,
                       stepFactor=1.5, 
                       improve=1e-5, ntree=15001)
      bestmtry.all <- rbind(bestmtry.all, bestmtry)
    }
     
     mtr <- bestmtry.all[bestmtry.all[,2] == min(bestmtry.all[,2]), ]
     
     if (class(mtr) == "numeric") {
       mtr.f <- mtr[1]
     } else {
         mtr.f <- mtr[1,1]
     }
     
     rf.m <- randomForest(Status ~., 
                                  data = rfdata.samp, 
                                  importance=TRUE,
                                  proximity=TRUE, 
                                  mtry=mtr,
                                  ntree=15001) 
    rf.models[[i]] <- rf.m
    
    rf.print <- capture.output(print(rf.m))
    
    samp.sammary[i, "OOBE_indRF"] <- as.character(sub(".*: ", "", rf.print[8]))
    
    samp.sammary[i, "Low_indRF"] <- rf.m$confusion[2,3]
    
    samp.sammary[i, "High_indRF"] <- rf.m$confusion[1,3]
    
    samp.sammary[i, "Mtry_indRF"] <- mtr.f
    
    samp.sammary[i, "Class_indRF"] <- as.character(predict(rf.m, samp))
    
  }




rfdataII.lean <- rfdataII[, colnames(rfdataII) %in% rf.lean.taxa ]

add.rfdataII.lean <- data.frame(matrix(data = 0, nrow = 94, ncol = 8))

colnames(add.rfdataII.lean) <- rf.lean.taxa[!rf.lean.taxa %in% colnames(rfdataII.lean)]

rfdataII.lean.final <- cbind(rfdataII.lean, add.rfdataII.lean)

pred.lean <- c()
for(i  in 1:94) {
  p <- predict(rf.lean.model, rfdataII.lean.final[i, ])
  pred.lean <- c(pred.lean, as.character(p))
  
}


sn <- data.frame(str_split(ps0.prop@sam_data$Sample.name, "_", simplify = TRUE)) 

sn$predition <- pred.lean


load("output/metadata.R")
shedding.data.t1.pl <- data.frame(t(shedding.data.t1))
shedding.data.t1.pl$Prediction <- pred.lean

rfdataII <- data.frame(ps0.prop@otu_table@.Data)

rf.full.taxa <- rownames(rf.full.model$importance)

rfdataII.full <- rfdataII[, colnames(rfdataII) %in% rf.full.taxa ]

miss.t.full <- rf.full.taxa[!rf.full.taxa %in% colnames(rfdataII.full)]

add.rfdataII.full <- data.frame(matrix(data = 0, nrow = 94, ncol = length(miss.t.full)))

colnames(add.rfdataII.full) <- rf.full.taxa[!rf.full.taxa %in% colnames(rfdataII.full)]

rfdataII.full.final <- cbind(rfdataII.full, add.rfdataII.full)

pred.full <- c()
for(i  in 1:94) {
  p <- predict(rf.full.model, rfdataII.full.final[i, ])
  pred.full <- c(pred.full, as.character(p))
  
}

sn$predition.full <- pred.full



samp.sammary$Class_lean <- sn$predition
samp.sammary$Class_full <- sn$predition.full



load("resource/rfdata.Rdata")

rfdata.seqnames <- rfdata

colnames(rfdata.seqnames) <- c(rf.full.taxa, "Status") 

sp1.II <- rfdataII[1,][,!colSums(rfdataII[1,]) == 0]

sp1.II.over <- colnames(rfdata.seqnames)[colnames(rfdata.seqnames) %in% colnames(sp1.II)]

rfdata.seqnames.sp1 <- rfdata.seqnames[ ,c(sp1.II.over, "Status")]

set.seed(395726)
rf.model.sp1 <- randomForest(Status ~., 
                              data = rfdata.seqnames.sp1, 
                              importance=TRUE,
                              proximity=TRUE, 
                              mtry=30,
                              ntree=15001) 


shedding.data.t2.m<- shedding.data.t1.m[shedding.data.t1.m$value %in% c("0", "+", "++", "+++"), ]

shedding.data.t2.m$value <- factor(shedding.data.t2.m$value, levels = c("0", "+", "++", "+++"))

shedding.data.t2.m$Score <- as.character(shedding.data.t2.m$value)

shedding.data.t2.m$Score[shedding.data.t2.m$Score %in% c("+", "++", "+++")] <- c(1,2,3)[match(shedding.data.t2.m$Score, c("+", "++", "+++"), nomatch = 0)]

samp.score <- aggregate(as.numeric(shedding.data.t2.m$Score), 
          by=list(Category=shedding.data.t2.m$Var2), FUN=sum)

samp.score <- samp.score[order(samp.score$x, decreasing = TRUE), ]

shedding.data.t2.m$Var2 <- factor(shedding.data.t2.m$Var2, levels = samp.score$Category)

shedding.plot <- ggplot(shedding.data.t2.m,
        aes(x=Var1, y=value, color=value)) + 
        geom_point(size=2) + 
        facet_wrap(~Var2, ncol = 3) + 
        theme_bw()  +               
        theme(axis.text.x = element_blank()) 

ggsave("output/test_plots/shedding_data.png", shedding.plot, dpi=300, width = 14, height = 8, units = "in")                              

samp.sammary$ID <- ps0.II.prop@sam_data$Sample.name
samp.sammary$Cow_ID <- sn$X4

samp.score$Category <- gsub("X", "", samp.score$Category)
rownames(samp.score) <- gsub("X", "", samp.score$Category)

for (sp in samp.score$Category) {
  samp.sammary$Sheddiing_score[sp == samp.sammary$Cow_ID] <- samp.score[sp,"x"]
}

samp.sammary <- samp.sammary[order(samp.sammary$Sheddiing_score, decreasing = TRUE), ]

write.csv(file = "output/test_plots/sam_tab.csv", samp.sammary)
