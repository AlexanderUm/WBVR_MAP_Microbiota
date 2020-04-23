#################################################
#Clean up metadata for MAP & Microbiota Run II##
################################################

#Set enviormnet
###############
setwd("C://Rprojects/Microbiota_MAP_II_20200417/data/metadata/")

set.seed(497495632)

library(reshape2)
library(tidyverse)


#Read in files
##############
shedding.data <- read.csv("all_shading_data_working_copy.csv")

seq.plate.data <- read.csv("20200323_v2_WBVR_BaseClearOrder_124511.csv")

#Format files 
#############
shedding.data.t1 <- shedding.data[,c(4:ncol(shedding.data))]

shedding.data.t1 <- apply(shedding.data.t1, 2, as.character)

shedding.data.t1[shedding.data.t1 %in% ""] <- "ND"

rownames(shedding.data.t1) <- gsub("-", "_", as.character(shedding.data$SDATE))

shedding.data.t1.m <- melt(shedding.data.t1)

shedding.data.t1.m$Sample.name <- paste0(shedding.data.t1.m$Var1, 
                                         gsub("X", "_", shedding.data.t1.m$Var2))

meta.clean.bII <- left_join(seq.plate.data, shedding.data.t1.m[,c(3,4)], by="Sample.name")

#Save objects
#############
save(list = c("shedding.data.t1.m", "meta.clean.bII"), 
          file = "C://Rprojects/Microbiota_MAP_II_20200417/output/metadata.R")
