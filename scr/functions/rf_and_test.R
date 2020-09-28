#Fucntion: build a RF model and test prediciton 
#################################################
# samples_to_test - will be substracted from Rf_data
# The model will be build once for rf_data without samples_to_test
# Prediction will be done n times (n_times_prediction) 

rf_and_test <- function (rf_data, samples_to_test, n_samples_training, mtry, ntree, variable_column, regression_TorF) {
    
    rf_data_f1 <- rf_data[!rownames(rf_data) %in% samples_to_test, ]
    
    test_samp <- rf_data[rownames(rf_data) %in% samples_to_test, ]
    
    if (regression_TorF == TRUE) {
        
        rf_data_f1[,variable_column] <- as.numeric(as.character(rf_data_f1[,variable_column]))
    }
    
    rf_model <- randomForest(as.formula(paste0(variable_column, " ~", " .")), 
                           sampsize = as.numeric(n_samples_training),  
                           data = rf_data_f1,
                           mtry=as.numeric(mtry),
                           ntree=as.numeric(ntree))
    
        
        pred_res <- predict(rf_model, test_samp)
        
    return(pred_res)
}