# Run RF model on defined proportion of randomly subseted dataset
#################################################################

rf_by_size <- function(data, sub_percent, class_colum, mtr, ntree, nrep, class_colors) {
    
# Prepare subsets of data for rf
    rf.data.l <- list()
    
    
    for (ssub in  sub_percent) { 
    
        for (r in 1:nrep) {
          
            rfdata.ssub <- c()
          
            for (gr in unique(data[, class_colum])) {
              
                data.sub1 <- data[data[, class_colum] %in% gr, ]
              
                data.sub2 <- data.sub1[sample(nrow(data.sub1), 
                                    round(ssub*nrow(data.sub1), 0), 
                                    replace = FALSE), ]
              
                 rfdata.ssub <- rbind(rfdata.ssub , data.sub2)
            }
            
           rfdata.ssub <- data.frame(rfdata.ssub) 
            
           df.name <- paste(ssub, "Alt", r, sep = "_")
            
           rf.data.l[[df.name]] <- rfdata.ssub
        }
    }
    
    # Build RF model for each subset  
    rf.list <- foreach (i=names(rf.data.l)) %dopar% {
        
            rf.obj <- randomForest::randomForest(as.formula(paste0(class_colum, " ~ .")), 
                                                data = rf.data.l[[i]], 
                                                importance=TRUE,
                                                proximity=TRUE, 
                                                mtry=mtr,
                                                ntree=ntree) 
        
            rf.obj$ID <- i
        
            rf.obj
            }
    
    # Extract data for plotting 
     ssub.data <- c()
    
     for (rf1 in 1:length(rf.list)) {
         
         rf <- rf.list[[rf1]]
         
         rf.out <- data.frame(rf$confusion)
         
         rf.out$ID <- rf$ID
         
         rf.out$Shedder <- rownames(rf.out)
         
         rf.out$Proportion <- as.numeric(gsub("_Alt_.*", "", rf$ID))
         
         ssub.data <- rbind(ssub.data, rf.out)
    }
    
    
    # Plot the data 
    ret.plot  <- ggplot(ssub.data, aes(y= class.error, x=Proportion, color = Shedder)) + 
              geom_jitter(data = ssub.data, aes(y= class.error, x=Proportion, color = Shedder), 
              width = 0.0025, height = 0.0025, size = 2, alpha=0.75, stroke = 1) + 
              stat_smooth(method = "loess", formula = y ~ x, size = 1) + 
              theme_bw() + 
              ylab("Error rate") + 
              xlab("Proportion of the data") + 
              scale_color_manual(values = class_colors)
              
    return(ret.plot)
    
    rm(list=ls())
    gc()
}
