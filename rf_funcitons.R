
###########################################################
#I. Find and plot best number of trees and mtry parameter #
###########################################################
###########################################################

TreeMtryPlot <- function (data, ntrees, class_colum, ntimes) {

  # Run the adjustment of mtry parrametar with various number of trees 
  #####################################################################
  mtry.out2<- list()

  for (tr in as.character(ntrees)) {
    mtry.out1 <- foreach(i=1:ntimes) %dopar% {randomForest::tuneRF(data[, -which(names(data) %in% class_colum)],
                                                                     data[, which(names(data) %in% class_colum)],
                                                                     round(sqrt(ncol(data)), 0), 
                                                                     ntreeTry=as.numeric(tr), 
                                                                     stepFactor=2, 
                                                                     improve=0.05,
                                                                     trace=TRUE, plot=TRUE, doBest=FALSE)
    }
    mtry.out2[[tr]] <-  mtry.out1
  }
  
  
  
  #Extract data into dataframe 
  #############################
  mtry.out2.df <- c()
  
  for (tr in as.character(ntrees)) {
    
    for (i in 1:ntimes) {
      
      one.try <- data.frame(mtry.out2[[tr]][[i]])
      
      one.try$ntrees <- tr
      
      one.try$run <- i
      
      mtry.out2.df <- rbind(mtry.out2.df, one.try)
    }
  }
  
  
  
  #Adjust levels in dataframe for cleaner plotting 
  ################################################
  mtry.out2.df$OOBError.round <- round(as.numeric(mtry.out2.df$OOBError), 2)
  
  mtry.out2.df$mtry <- factor(mtry.out2.df$mtry, 
                              levels = unique(mtry.out2.df$mtry)[order(unique(mtry.out2.df$mtry))])
  
  mtry.out2.df$ntrees <- factor(mtry.out2.df$ntrees, 
                                levels = unique(mtry.out2.df$ntrees)[order(
                                  as.numeric(
                                    as.character(unique(mtry.out2.df$ntrees))))])
  
  
  #Plot results
  #############
  plot.mtr.ntr <- ggplot(data = mtry.out2.df, aes(x=mtry, y=OOBError.round, group=run, color=factor(run))) + 
    geom_point(size=3, alpha = 0.5) + 
    geom_line(alpha = 0.7) +
    facet_grid(ntrees ~., scales = "free") + 
    theme_bw() 

  #Returnt the plot 
    return(plot.mtr.ntr) }
  

#####################################
#II. Accuracy per split             #
#####################################
#####################################


BestSplitAccuracy <- function(data, ntrees, class_colum, mtry, ntimes) {
  
     out1 <- list()
    
     for (m in as.character(mtry)) {
      
       out2 <- list()
       
        for (nt in as.character(ntrees)) {
          
          rf.rep <- foreach(i=1:ntimes) %dopar% {randomForest::randomForest(as.formula(paste0(class_colum, " ~ .")), 
                                                                 data = data.frame(data), 
                                                                 importance=TRUE,
                                                                 proximity=TRUE, 
                                                                 mtry=as.numeric(m),
                                                                 ntree=as.numeric(nt))
            }
          
          out2[[nt]] <- rf.rep
        }
        
     out1[[m]] <- out2
     
     }
     
      

     
#Extract information about class errors 
     
     rf.ntree <- c()
     
     for (rf1 in names(out1)) {
       
       out1.rf1 <- out1[[rf1]]
       
       for (rf2 in names(out1.rf1)) {
         
         out1.rf2 <- out1.rf1[[rf2]]
         
         for (rf3 in 1:ntimes) {
           
            rf.out <- data.frame(out1.rf2[[rf3]]$confusion)
            
            rf.out$mtry <- rf1
            
            rf.out$ntree <- rf2
            
            rf.out$atempt <- rf3
            
            rf.out$Status <- rownames(rf.out)
            
            rf.ntree <- rbind(rf.ntree, rf.out)
         
            } 
         
       }
       
    }

#Convert informaiton about class errors into dataframe 
rf.ntree.df <- data.frame(rf.ntree)

rf.ntree.df$ntree <- factor(rf.ntree.df$ntree, 
                            levels = unique(rf.ntree.df$ntree)[order(
                              as.numeric(
                                as.character(unique(rf.ntree.df$ntree))))])


#Plot class errors 
best.tree.split <- ggplot(rf.ntree.df, aes(y= as.numeric(as.character(class.error)),
                                           x=ntree, color = Status)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1, height = 0.001, size = 3, alpha=0.75) + 
  theme_bw() + 
  facet_grid(~ mtry) +
  ylab("Error rate")

return(best.tree.split)

}



###################################################################
#III. Influence of the dataset size on accuracy of classification #
###################################################################
###################################################################

#3.1 Define function that will subest defined presetage of samples from High and Low (each group)
#    shading groups 





PlotSizeInfluece <- function(data, sub.percent, class_colum, mtr, ntree, nrep) {
  
  
err.l <- list()

for (ssub in  sub.percent) { 
  
  rf.rep <- foreach(i=1:nrep) %dopar% {
 
       
 # Prepare data for rf 
      rfdata.ssub <- c()
    
      for (gr in levels(data[, class_colum])) {
      
        data.sub1 <- data[data[, class_colum] %in% gr, ]
      
        data.sub2 <- data.sub1[sample(nrow(data.sub1), 
                                    round(ssub*nrow(data.sub1), 0), 
                                    replace = FALSE), ]
      
        rfdata.ssub <- rbind(rfdata.ssub , data.sub2)
                                              }
    
    
 # Random forest function 
    randomForest::randomForest(as.formula(paste0(class_colum, " ~ .")), 
                                                                data = data.frame(rfdata.ssub), 
                                                                importance=TRUE,
                                                                proximity=TRUE, 
                                                                mtry=mtr,
                                                                ntree=ntree) 
    
        } 


err.l[[as.character(ssub)]] <- rf.rep

}
  


# Extract data from the resulted object 
ssub.data <- c()
for (rf1 in names(err.l)) {
  for (rf2 in 1:nrep) {
    rf.out <- data.frame(err.l[[rf1]][[rf2]]$confusion)
    rf.out$ntree <- rf1
    rf.out$Status <- rownames(rf.out)
    
    ssub.data <- rbind(ssub.data, rf.out)
  }
}

# Plot class errors per group with fitted line (Locally Weighted Regression)
ssub.plot <- ggplot(ssub.data, aes(y= as.numeric(as.character(class.error)),
                                   x=as.numeric(ntree), color = Status)) + 
  geom_jitter(width = 0.01, height = 0.001, size = 2, alpha=0.75) + 
  stat_smooth(method = "loess", formula = y ~ x, size = 1) + 
  theme_bw() + 
  ylab("Error rate") + 
  xlab("Proportion of the data")

return(ssub.plot)

}


  