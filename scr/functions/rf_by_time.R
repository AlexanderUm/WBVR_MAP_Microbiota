# Functinon that will take data and slice dataset by order of variables in suplied metadata. 
# #########################################################################################

rf_by_time <- function(data_for_rf, metadata_for_rf, class_column, slice_column, 
                       slice_steps, mtr, ntree, nrep, class_colors) {
    
    # Prepare data
    sliced.datasets <- list()
    
    data.for.slice <- data_for_rf
    
    metadata.for.slice <- metadata_for_rf
    
    # Slice data by supplied array 
    for (slice in slice_steps) {
        
        metadata.for.slice <- metadata.for.slice[!metadata.for.slice[,slice_column] %in% slice, ]
        
        data.for.slice <- data.for.slice[rownames(metadata.for.slice), ]
        
        sliced.datasets[[slice]] <- data.for.slice }

    # Replicate n times 
    slice_steps.rep <- rep(sliced.datasets, nrep) 
    
    # Run RF models for each "sliced" data set 
    sliced.rf <- foreach (s=1:length(slice_steps.rep), .packages = "randomForest") %dopar% { 
        
                            rf.obj <- randomForest(as.formula(paste0(class_column, " ~ .")), 
                                      data = slice_steps.rep[[s]], 
                                      importance=FALSE,
                                      proximity=FALSE, 
                                      mtry=mtr,
                                      ntree=ntree)
        
                            rf.obj$Sliced <- names(slice_steps.rep)[s]
        
                            rf.obj$ncol <- nrow(slice_steps.rep[[s]])   
        
                            rf.obj
                            }
    
    # Extract data for plotting 
    ssub.data <- c()
    
    for (rf1 in 1:length(sliced.rf)) {
                        
        rf <- sliced.rf[[rf1]]
                        
        rf.out <- data.frame(rf$confusion)
                        
        rf.out$Slice <- rf$Slice
                        
        rf.out$Shedder <- rownames(rf.out)
                        
        rf.out$Proportion <- round(length(rf$y)/nrow(metadata_for_rf), 2)
        
        ssub.data <- rbind(ssub.data, rf.out)
        }
    
    ssub.data$Slice <- factor(ssub.data$Slice, levels=c(slice_steps))   
              
    # Plot 
    ret.plot  <- ggplot(ssub.data, aes(y= class.error, x=Proportion, color = Shedder)) + 
              geom_jitter(data = ssub.data, aes(y= class.error, x=Proportion, color = Shedder, shape=Slice), 
              width = 0.0025, height = 0.0025, size = 2, alpha=0.75, stroke = 1.5) + 
              stat_smooth(method = "loess", formula = y ~ x, size = 1) + 
              theme_bw() + 
              scale_shape_manual(values=c(1:length(unique(ssub.data$Slice)))) + 
              ylab("Error rate") + 
              xlab("Proportion of the data")  + 
              scale_color_manual(values = class_colors)
             
    return(ret.plot)
    
    rm(list=ls())
    gc()
}
