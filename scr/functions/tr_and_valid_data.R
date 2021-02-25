# Function to prepare training and validataion datasets for "early", 
# "middle", and "late" groups

tr_and_valid_data <- function (samp_valid, tr_subset_prop, rf_data_in) {

    data.set.size <- floor(length(samp_valid)/tr_subset_prop)

    # Generate a random sample of "data.set.size" indexes
    indexes <- samp_valid[sample(1:length(samp_valid), size = data.set.size)]

     # Assign the data to the correct sets
    training <- rf_data_in[!rownames(rf_data_in) %in% indexes,]

    validation <- rf_data_in[ indexes,]
    
    out <- list(training = training, 
                validation = validation)
    
    return(out)
   
}