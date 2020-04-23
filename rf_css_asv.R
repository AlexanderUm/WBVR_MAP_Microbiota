
#0. Prepare enviorment 
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

#Prepare data 

css.rf.data <- data.frame(t(ps1.1.css@otu_table))
css.rf.data$Status <- ps1.1.css@sam_data$Status


################################################
# 2. Optimization of parameters for RF model 
################################################

#2.1 Optimiztion of nubber of features used per split used for with different number of trees. 
##################################################################################

#2.1.1 list number of used trees 
ntrees <- c(1000, 5000, 10000, 15000, 20000)

#2.1.2 Create an empty list 
mtry.out2.css <- list()

#2.1.3 Run the adjustment of mtry parrametar with various number of trees 
#Set seed 
set.seed(4574)

for (tr in as.character(ntrees)) {
  mtry.out1 <- foreach(i=1:7) %dopar% {  randomForest::tuneRF(css.rf.data[, -which(names(css.rf.data) %in% "Status")],
                                                              css.rf.data[, which(names(css.rf.data) %in% "Status")],
                                                              44, ntreeTry=as.numeric(tr), stepFactor=2, improve=0.05,
                                                              trace=TRUE, plot=TRUE, doBest=FALSE)
  }
  mtry.out2[[tr]] <-  mtry.out1
}
rm(list = c("tr", "i", "mtry.out1"))

#mtry.out2.css <- mtry.out2
#2.1.4 Save resulted object 
#save(mtry.out2.css, file="output/rf_objects/mtry_out2_css.RData")

#2.1.5 Convert the object in to dataframe 
mtry.out2.df <- c()
for (tr in as.character(ntrees)) {
  for (i in 1:7) {
    one.try <- data.frame(mtry.out2[[tr]][[i]])
    one.try$ntrees <- tr
    one.try$run <- i
    
    mtry.out2.df <- rbind(mtry.out2.df, one.try)
  }
}
rm(list = c("tr", "i", "one.try", "mtry.out2"))
gc()

#2.1.6 Adjust dataframe for correct plotting 
mtry.out2.df$OOBError.round <- round(as.numeric(mtry.out2.df$OOBError), 2)

mtry.out2.df$mtry <- factor(mtry.out2.df$mtry, 
                            levels = c("11", "22", "44", "88", "176", "352", "704"))
mtry.out2.df$ntrees <- factor(mtry.out2.df$ntrees, 
                              levels = c("1000", "5000", "10000", "15000", "20000"))


#2.1.7 Plot results 
plot.mtr.ntr <- ggplot(data = mtry.out2.df, aes(x=mtry, y=OOBError.round, group=run, color=run)) + 
  geom_point(size=3, alpha = 0.5) + 
  geom_line(alpha = 0.7) +
  facet_grid(ntrees ~., scales = "free") + 
  theme_bw()

#2.1.8 Save the plot 
ggsave("output/processing/ntree_mtry_asv.pdf", plot = plot.mtr.ntr)
rm(list = c("plot.mtr.ntr", "mtry.out2.df"))


#2.2 Buld RF models in 8 folds using 176 and 352 fetures per split (best parameters)

#2.2.1 Buld 8 RF models with each number of trees and mtry=176
ntrees.out.176 <-list()
for(nt in as.character(ntrees)) {
  rf.rep <- foreach(i=1:8) %dopar% {  randomForest::randomForest(Status ~ ., 
                                                                 data = css.rf.data, 
                                                                 importance=TRUE,
                                                                 proximity=TRUE, 
                                                                 mtry=176,
                                                                 ntree=as.numeric(nt))
  }
  ntrees.out.176[[nt]] <- rf.rep
  gc()
}

#2.2.2 Buld 8 RF models with each number of trees and mtry=352
ntrees.out.352 <-list()
for(nt in as.character(ntrees)) {
  rf.rep <- foreach(i=1:8) %dopar% {  randomForest::randomForest(Status ~ ., 
                                                                 data = css.rf.data, 
                                                                 importance=TRUE,
                                                                 proximity=TRUE, 
                                                                 mtry=352,
                                                                 ntree=as.numeric(nt))
  }
  ntrees.out.352[[nt]] <- rf.rep
  gc()
}

#2.2.3 Remove temp loop objects 
rm(list=c("nt", "ntrees", "rf.rep"))

#2.2.4 Save resulted object 
save(list = c("ntrees.out.176", 
              "ntrees.out.352"), 
     file = "output/rf_objects/ntrees.out.RData")


#2.2.5 Extract information about class errors 
rf.ntree <- c()
for (rf1 in names(ntrees.out.176)) {
  for (rf2 in 1:8) {
    rf.out <- data.frame(ntrees.out.176[[rf1]][[rf2]]$confusion)
    rf.out$ntree <- rf1
    rf.out$nmtry <- "splist@176"
    rf.out$Status <- rownames(rf.out)
    
    rf.out2 <- data.frame(ntrees.out.352[[rf1]][[rf2]]$confusion)
    rf.out2$ntree <- rf1
    rf.out2$nmtry <- "splist@352"
    rf.out2$Status <- rownames(rf.out2)
    
    rf.ntree <- rbind(rf.ntree, rf.out, rf.out2)
  }
}
rm(list=c("rf1", "rf1.out", "rf2", "rf.out2", "ntrees.out.176", "ntrees.out.352"))
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
ggsave("output/processing/best_tree_split.pdf", plot = best.tree.split)                       
rm(list=c("rf1", "rf1.out", "rf2", "rf.out2")) 
