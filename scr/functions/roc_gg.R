# Function to plot ROC and culculate ROC
########################################
roc_gg <- function(train_set, valid_set, mtry, ntree, class_col) {
    
    # Build rf model 
    lean.rf <- randomForest(as.formula(paste0(class_col, " ~ .")), 
             mtry = mtry,
             data = train_set, 
             importance = TRUE,
             proximity = TRUE, 
             ntree = ntree) 
    
    # Predict values for validation dataset 
    predictions <- as.numeric(predict(lean.rf, valid_set[, !colnames(valid_set) %in% class_col]))
    
    # Build ROC                                                     
    pred <- prediction(predictions, valid_set[, colnames(valid_set) %in% class_col])
    
    perf <- performance(pred, measure = "tpr", x.measure = "fpr")
     
    # Culculate AUC                                          
    auc.perf <- performance(pred, measure = "auc")
                                              
    # Combine data to report it                                           
    data.out <- data.frame(cbind(perf@y.values[[1]], perf@x.values[[1]], auc.perf@y.values[[1]]))

    colnames(data.out) <- c("y_values", "x_values", "AUC")

    oob <- capture.output(lean.rf)

    oob <- oob[grep("OOB" ,oob)]

    data.out$OOBE <- trimws(sub(".*:", "", oob))

    class.err <- lean.rf$confusion[,"class.error"]

    for (i in 1:length(class.err)) {
        data.out[,paste0(names(class.err[i]), "_error")] <- class.err[i]
        }
    
    return(data.out)
   
}