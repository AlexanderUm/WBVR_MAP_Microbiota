###########################################################
#I. Find and plot best number of trees and mtry parameter #
###########################################################
###########################################################

Tree_Mtry_Plot <- function (data, ntrees, start_val, class_colum, ntimes, stepF) {

  # Run the adjustment of mtry parrametar with various number of trees 
  #####################################################################
  mtry.out2<- list()

  for (tr in as.character(ntrees)) {
    mtry.out1 <- foreach(i=1:ntimes) %dopar% {randomForest::tuneRF(data[, -which(names(data) %in% class_colum)],
                                                                     data[, which(names(data) %in% class_colum)],
                                                                     mtryStart = start_val, 
                                                                     ntreeTry=as.numeric(tr), 
                                                                     stepFactor=stepF, 
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