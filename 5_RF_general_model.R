##############################################################
# General random forest model 
# - Find the best assign shedding class based on accuracy 
#   of classification. 
# - Optimize the model parameters 
##############################################################

# Set seed
set.seed(4597234)

# Load required libraries. 
source("scr/functions/general/load_abs_install_pkg.R")

load_abs_install_pkg(c("phyloseq", "tidyverse", "foreach", "doParallel",
                       "randomForest", "rfPermute", "ComplexHeatmap"))

# Create output directory. 
dir.create("output/5_RF_general_model")


# 1. Extract and format data for RF model building 
##################################################

# Read in filtered and normalized phyloseq object 
ps.tf.css.01 <- readRDS("output/3_filt_norm_phyl/ps_tf2_css.RDS")

# Load custom function to format data from phyloseq for RF analysis 
source("scr/functions/data_for_rf.R")

# Format data for RF regression analysis
rf.data.css.01 <- data_for_rf(phyloseq = ps.tf.css.01, 
                              class.column = 'WeightedScoreII', 
                              remove.taxa.prev.less.than = 1, 
                              return.df = TRUE)

# Convert the WeightedScoreII column into numeric 
rf.data.css.01$WeightedScoreII <- as.numeric(as.character(rf.data.css.01$WeightedScoreII))



# 2. Create a list of datasets with different configuration of shedding groups
##############################################################################

# Extract shedding score from data for RF
score <- unique(rf.data.css.01$WeightedScoreII)

# Order shedding scores
score <- score[order(score, decreasing = TRUE)]

# Make an empty list for datasets with different configuration of shedding groups
rfd.cat.l <- list()

# Create datasets with different configuration of High and Low groups (loop)
# "s" is a unique weighted shedding score
# First and last 4 shedding scores are skipped to prevent dramatic differences
#           between "Low" and "High" shedding groups
for (s in 4:(length(score)-4)) {
    
    # Copy RF data object
    rfd.cat <- rf.data.css.01
    
    # Change the name of the response column
    colnames(rfd.cat)[colnames(rfd.cat) %in% "WeightedScoreII"] <- "Shedder"
    
    # Assign samples to "Low" or "High" shedding group base on the condition
    rfd.cat$Shedder <- as.factor(ifelse(rfd.cat$Shedder < score[s], "Low", "High"))
    
    # Add add RF data to list
    rfd.cat.l[[s]] <- rfd.cat
}


# 3. Build RF models 
####################

# Make a computing cluster (number of cores)
cl <- makeCluster(10)

# Register the cluster
registerDoParallel(cl)

# Run RF models building in parallel 
rf.cat.l.res <- foreach(p = 4:length(rfd.cat.l),.packages = "randomForest") %dopar% {
    
                randomForest(Shedder ~ ., 
                             data=rfd.cat.l[[p]], 
                             importance=TRUE,
                             proximity=TRUE, 
                             ntree=7501) 
                }

# Stop the computing cluster 
stopCluster(cl)


# 4. Extract and format data for visualization 
##############################################

# Create empty vector for formatted data
rf.plot.d <- c()

# Extraction of RF confusion data from RF objects (loop)
# "i" is number of RF object in the list 
for (i in 1:length(rf.cat.l.res)) {
    
    # Extract RF object from the list 
    rf.ind <- rf.cat.l.res[[i]]
    
    # Combine data from RF object confusion matrix into a long matrix    
    rf.plot.d <- rbind(rf.plot.d, cbind(rf.ind$confusion, rep(i, 2)))
}

# Covert combine confusion matrix into dataframe 
rf.plot.d <- data.frame(rf.plot.d)

# Remove dots "." from row names and create from it column "Class"
rf.plot.d$Class <- sub("\\.", "", rownames(rf.plot.d))

# Remove numbers from the class column 
rf.plot.d$Class <- gsub("[0-9]", "", rf.plot.d$Class)

# Add column "Split", that reflects weighted score that was used for split 
#             between "Low" and "High" shedders 
rf.plot.d$Split <- as.character(rep(score[4:(length(score)-4)], each = 2))

# Write the combine confusion table into a file. This file will be use 
#           for plotting of results in the Section 7 
write.csv(rf.plot.d, "output/5_RF_general_model/best_split.csv")


# 5. Plot class errors (confusion matrix) at different composition of shedding groups
####################################################################################

split.rf <- ggplot(rf.plot.d, aes(x=Split, y=class.error, color=Class)) + 
                geom_point(size=3) + 
                theme_bw() + 
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +                
                scale_color_manual(values = c("steelblue", "gold3"))

# Save the plot 
ggsave(filename = "output/5_RF_general_model/Figure_3_1.png", 
       plot = split.rf, width = 5, height = 3.5, dpi = 400)
ggsave(filename = "output/5_RF_general_model/Figure_3_1.pdf", 
       plot = split.rf, width = 5, height = 3.5)


# 6. Determine optimal parameters for RF model with full dataset  
################################################################

# Copy "rf.data.css.01" object into a new object 
rf.data.css <- rf.data.css.01

# Change name of the "WeightedScoreII" to "Shedder"
colnames(rf.data.css)[colnames(rf.data.css) %in% "WeightedScoreII"] <- "Shedder"

# Define all samples with Shedding score less than 0.51 as "Low" Shedders 
rf.data.css$Shedder <- as.factor(ifelse(rf.data.css$Shedder < 0.51, "Low", "High"))

# Load custom function for Mtry tuning and visualization 
source("scr/functions/Tree_Mtry_Plot.R")

# Make a computing cluster (number of cores)
cl <- makeCluster(10)

# Register the cluster 
registerDoParallel(cl)

# Visualize Mtry tunning attempts 
tree.mtry.plot.all <- Tree_Mtry_Plot(data = rf.data.css, 
                                     ntrees = c(7001, 10001, 15001), 
                                     start_val = round(sqrt(ncol(rf.data.css)), 0), 
                                     stepF = 0.5, 
                                     class_colum = "Shedder", 
                                     ntimes = 5)

# Stop the computing cluster
stopCluster(cl)

# Save the plot 
ggsave(filename = "output/5_RF_general_model/Figure_S6.png", plot = tree.mtry.plot.all, 
       width = 5, height = 4,  dpi = 400) 
ggsave(filename = "output/5_RF_general_model/Figure_S6.pdf", plot = tree.mtry.plot.all, 
       width = 5, height = 4) 
