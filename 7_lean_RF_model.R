#################
# "Lean" RF model
#################

# Set environment
#################

# Set seed 
set.seed(49562)

# Load function for loading and installing libraries 
source("scr/functions/general/load_abs_install_pkg.R")

# Load libraries 
load_abs_install_pkg(c("ggpubr", "phyloseq", "tidyverse", "foreach", 
                       "doParallel", "randomForest", "ROCR", "RColorBrewer", 
                       "gmodels"))

# Create directory 
dir.create("output/7_lean_RF_model")


# 1. Prepare data for the "Lean" RF model 
#########################################

# Extract and format data for RF from phyloseq using a custom function 
# Load phyloseq object
ps.tf.css.01 <- readRDS("output/3_filt_norm_phyl/ps_tf2_css.RDS")

#  Load custom function for formatting data for RF
source("scr/functions/data_for_rf2.R")

# Format the data 
rf.data.css.01 <- data_for_rf2(phyloseq = ps.tf.css.01, 
                              class.column = 'WeightedScoreII', 
                              remove.taxa.prev.less.than = 1)

# Subset only significantly contributing taxa 
# Read table with significantly contributed taxa from generated in script 6. 
sig.tax <- read.table("output/6_discriminatory_taxa/sig_contr_taxa.txt")

# Subset only significantly contributing taxa from the data for RF
rf.data.css.sig <- rf.data.css.01[, rownames(sig.tax)]

# Add Weighted Score column 
rf.data.css.sig$WeightedScoreII <- rf.data.css.01$WeightedScoreII

# Create a list of datasets with different configuration of the shedding groups 
# Make a vector of unique shedding scores 
score <- unique(rf.data.css.sig$WeightedScoreII)

# Order unique shedding scores in decreasing order 
score <- score[order(score, decreasing = TRUE)]

#  Make an empty list for variously configured RF data 
rfd.cat.l <- list()

# "s" stands for an entry in the vector with unique shedding scores.
# first and last 4 shedding scores will not be taken into account.
for (s in 4:(length(score)-4)) {
    
    # Clon object 
    rfd.cat <- rf.data.css.sig
    
    # Rename column WeightedScoreII into Shedder 
    colnames(rfd.cat)[colnames(rfd.cat) %in% "WeightedScoreII"] <- "Shedder"
    
    # If shedding core less than "s" assign sample as "Low" shedder 
    rfd.cat$Shedder <- as.factor(ifelse(rfd.cat$Shedder < score[s], "Low", "High"))
    
    # Add dataset to list of the datasets 
    rfd.cat.l[[s]] <- rfd.cat
    
}

# 2. Make serious of "Lean" RF models 
#####################################

# Register a cluster for the multicore computing 
cl <- makeCluster(6)
registerDoParallel(cl)

# Run in parallel 
rf.cat.l.res <- foreach(p=4:length(rfd.cat.l),.packages = "randomForest") %dopar% {
    
                randomForest(Shedder ~ ., 
                             data=rfd.cat.l[[p]], 
                             importance=TRUE,
                             proximity=TRUE, 
                             ntree=7501) 
                }

# Stop cluster 
stopCluster(cl)


# 3. Visulaze class error per created "Lean" RF model 
#              with various configurations of Shedding groups
#############################################################

# Prepare data for plotting 
# Make empty object (vector) to extract data from the list of RF model into in 
rf.plot.d <- c()

# "i" RF object number in the list
for (i in 1:length(rf.cat.l.res)) {
    
    # Extract an individual RF object 
    rf.ind <- rf.cat.l.res[[i]]
    
    # Bind confusion matrices into a long table
    rf.plot.d <- rbind(rf.plot.d, cbind(rf.ind$confusion, rep(i, 2)))
}

# Convert into dataframe 
rf.plot.d <- as.data.frame(rf.plot.d)

# Format and add Class column 
 rf.plot.d$Class <- sub("\\.", "", rownames(rf.plot.d))

# Format Class column 
rf.plot.d$Class <- gsub("[0-9]", "", rf.plot.d$Class)

# Add Split column (score that was used for split into "High" and "Low" shedders)
 rf.plot.d$Split <- as.character(rep(score[4:(length(score)-4)], each = 2))

# Plot and class error per at each split category  

lean.rf <- ggplot(rf.plot.d, aes(x=Split, y=class.error, color=Class)) + 

                                         geom_point(size=3) + 

                                         theme_bw() + 

                                         theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + 
                                        
                                         scale_color_manual(values = c("steelblue", "gold3"))

ggsave(filename = "output/7_lean_RF_model/Figure_3_2.png", plot = lean.rf, width = 5, height = 3.5, dpi = 300)
ggsave(filename = "output/7_lean_RF_model/Figure_3_2.pdf", plot = lean.rf, width = 5, height = 3.5)


# 4. Combined visualization of class error "Full" and "Lean" RF model 
#              with various configurations of Shedding groups
######################################################################

# Load data from the general RF model 
split.d.gereal.RF <- read.csv("output/5_RF_general_model/best_split.csv")

# Add column Model_type 
split.d.gereal.RF$Model_type <- "General RF Model"

rf.plot.d$Model_type <- "Lean RF Model"

# Combine plotting data 
rf.comb.d <- rbind(split.d.gereal.RF[,-1], rf.plot.d)

# Format the split column  
rf.comb.d$Split <- gsub(",", ".", rf.comb.d$Split)

# Plot and save 
rf.split.comb <- ggplot(rf.comb.d, aes(x = Split, y = class.error, color = Class, group = Class)) + 
                                         geom_point(size = 2.5) + 
                                         geom_line() +
                                         theme_bw() + 
                                         theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
                                         scale_color_manual(values = c("steelblue", "gold3")) + 
                                         facet_grid( ~ Model_type) + 
                                         xlab("Split point (Weighted SIS)") + 
                                         ylab("Class Error")

ggsave(filename = "output/7_lean_RF_model/Figure_3.png", plot = rf.split.comb, width = 9, height = 3.5, dpi = 400)
ggsave(filename = "output/7_lean_RF_model/Figure_3.pdf", plot = rf.split.comb, width = 9, height = 3.5)

# 5. Optimization of "Lean" RF parameters (mtry and ntree)
##########################################################

# Prepare data
# Rename WeightedScoreII column into Shedder column 
colnames(rf.data.css.sig)[colnames(rf.data.css.sig) %in% "WeightedScoreII"] <- "Shedder" 

# Assign samples with Weighted Score (Shedding) less than 0.51 as "Low" shedders
rf.data.css.sig$Shedder <- as.factor(ifelse(rf.data.css.sig$Shedder < 0.51, "Low", "High"))

# Load function to test ntree and mtry (custom).
source("scr/functions/Tree_Mtry_Plot.R")

# Register a cluster that will be used by the function 
cl <- makeCluster(5)
registerDoParallel(cl)

# Run the function 
tree.mtry.plot.all <- Tree_Mtry_Plot(data = rf.data.css.sig, 
                                     ntrees = c(7501, 10001, 15001), 
                                     start_val = round(sqrt(ncol(rf.data.css.sig)), 0), 
                                     stepF = 0.5, 
                                     class_colum = "Shedder", 
                                     ntimes = 5)

# Stop the cluster 
stopCluster(cl)

# Save resulted plots 
ggsave(filename = "output/7_lean_RF_model/Figure_S8.png", plot = tree.mtry.plot.all, dpi = 400) 
ggsave(filename = "output/7_lean_RF_model/Figure_S8.pdf", plot = tree.mtry.plot.all) 


# 6. Test accuracy of lean RF model  
###################################

# Prepare data for calculation of prediction accuracy (AUC)
#         for complete and age groups datasets.

# Colmite dataset 
all.samp <- rownames(rf.data.css.sig)

# Age groups datasets
# Extract samples data from phyloseq to separate samples in groups by age
samp.meta <- sample_data(ps.tf.css.01)

# Early life (younger 12 month)
early.samp <- rownames(samp.meta[samp.meta$AgeMonth < 12, ])

# Middle life (12 to 24 month old)
middle.samp <- samp.meta[samp.meta$AgeMonth > 12, ]

middle.samp <- rownames(middle.samp[middle.samp$AgeMonth < 24, ])

# Late life (older than 24 month)
late.samp <- rownames(samp.meta[samp.meta$AgeMonth > 24, ])

# Combine datasets into a list 
roc.samp <- list(All = all.samp, 
                 Early = early.samp, 
                 Middle = middle.samp, 
                 Late = late.samp)

# Make a vector with denominators for calculation of number of samples
#        in training datasets.
prop.samp.roc <- c(5, 2, 2, 2)


# 7. Build RF and determine (ROC)AUC for each group 
###################################################

# Load custom function for running RF and data preparation for
#      visualisation of results as ROC.
source("scr/functions/roc_gg.R")

# Load custom function for preparation of training and validation data  
source("scr/functions/tr_and_valid_data.R")

# Make cluster 
cl <- makeCluster(10)
registerDoParallel(cl)

# Create empty vector to store data for visualisation 
roc.d.comb <- c()

# "d" object number in list of samples 
for (d in 1:length(roc.samp)) {
    
    # Start a parallel loop 
    roc.d.all <- foreach(i = 1:99, .packages = c("randomForest", "ROCR"))  %dopar%  { 
    
        # Prepare training and validation data with the custom function
        data.roc <- tr_and_valid_data(samp_valid = roc.samp[[d]],
                                 tr_subset_prop = prop.samp.roc[d],
                                 rf_data_in = rf.data.css.sig)
      
        # Prepare data for ROC visualisation
        roc_gg(train_set = data.roc$training, 
                   valid_set = data.roc$validation, 
                   mtry = 28, 
                   ntree = 15001, 
                   class_col = "Shedder")
        }
    
    # Add the Run column for each dataset with the Run number
    for (i in 1:length(roc.d.all)) {
    
        roc.d.all[[i]]$Run <- paste0("Run_", i)
        }
    
    # Combine into a long dataframe 
    roc.d.all <- do.call(rbind, roc.d.all)
    
    # Add column group with the name of the used dataset 
    roc.d.all$Group <- names(roc.samp)[d]
    
    # Bind into a single dataframe 
    roc.d.comb <- rbind(roc.d.comb, roc.d.all)
}

# Stop cluster 
stopCluster(cl)

# Write table into a file 
write.table(roc.d.comb, "output/7_lean_RF_model/roc_auc.txt")


# 8. Calculate CI and SD for AUC 
#################################

# Make an empty vector for confidence intervals 
ci.auc <- c()

# Calculate confidence interval per group 
# "g" group ID 
for (g in unique(roc.d.comb$Group)) {
    
    # Subset only samples belonging to the group 
    cd <- roc.d.comb$AUC[roc.d.comb$Group %in% g]
    
    # Calculate CI, SD and combine into a matrix 
    ci.auc <- rbind(ci.auc, c(ci(cd, confidence = 0.95), sd(cd), g))
    
}

# Change column names 
colnames(ci.auc) <- c("Estimate", 'CI lower', 'CI upper', 'Std. Error', "Std.Dev", "Group")

# Save as a table 
write.table(ci.auc, "output/7_lean_RF_model/auc_stat.txt")


# 9. Plot ROC for "Lean" RF 
############################

# Format the data for plotting
roc.d.comb$GroupUnique <- paste(roc.d.comb$Group, roc.d.comb$Run, sep = "_")

# Arrange factors for plotting in the Group column 
roc.d.comb$Group <- factor(roc.d.comb$Group, levels = c("All", "Early", "Middle", "Late"))

# Make a color palette for plotting
roc.d.colors <- brewer.pal(4, "Set1")

# Make an extra dataframe with median coordinates for ROC curves 
# Calculate means per group 
roc.d.means <- aggregate(roc.d.comb[, 1:2], list(roc.d.comb$Group), median)

# Make an empty object for the final dataframe 
roc.d.means.f <- c()

# "i" group name 
for (i in unique(roc.d.means$Group.1)) {
    
    # Bind mean value, [0,0], and [1,1] coordinates 
    mean.plus <- rbind(c(i, 0, 0), roc.d.means[roc.d.means$Group.1 %in% i, ], c(i, 1, 1))
    
    # Bind rows into a matrix 
    roc.d.means.f <- rbind(roc.d.means.f, mean.plus)  
}

# Convert x columns to numeric 
roc.d.means.f$x_values <- as.numeric(gsub(",", "\\.", roc.d.means.f$x_values))

# Convert y column to numeric 
roc.d.means.f$y_values <- as.numeric(gsub(",", "\\.", roc.d.means.f$y_values))

# Adjust column names 
colnames(roc.d.means.f) <- c("Group", "y_values", "x_values")

# Create an extra vector for colors 
roc.d.colors1 <- as.vector(t(replicate(3, roc.d.colors)))

# Plot and save ROC curve graph 
roc.all <- ggplot(roc.d.comb, aes(y = y_values, x = x_values, color = Group, group = GroupUnique)) + 
       geom_line(size = 0.2, alpha = 0.15) + 
       geom_abline(intercept = 0, slope = 1) + 
       theme_bw() + 
       ylab("True Positive rate") + 
       xlab("False Positive rate") +
       scale_color_manual(values = roc.d.colors, 
                          guide = guide_legend(override.aes = list(size = 1, alpha = 1)) ) + 
       geom_line(inherit.aes = FALSE, data = roc.d.means.f, 
                 mapping = aes(y = y_values, x = x_values, 
                               group = Group), color = roc.d.colors1, size = 0.7) 

ggsave(filename = "output/7_lean_RF_model/Figure_6B.png", plot = roc.all, width = 5, height = 3.5, dpi = 400)

ggsave(filename = "output/7_lean_RF_model/Figure_6B.pdf", plot = roc.all, width = 5, height = 3.5)


# 10. Make a AUC boxplot with significance levels
#################################################

# Remove redundant rows 
roc.d.box <- roc.d.comb[!duplicated(roc.d.comb$GroupUnique), ]

# Make a list of groups pears for comparison
my_comparisons <- list( c("All", "Early"), c("Early", "Middle"), c("Middle", "Late"), c("Early", "Late"), c("All", "Late"))

# Plot and save the box plot 
auc.box <- ggboxplot(roc.d.box, x = "Group", y = "AUC",
                color = "Group", palette = roc.d.colors,
                add = "jitter") + 
                stat_compare_means(comparisons = my_comparisons, label = "p.signif") 


ggsave(filename = "output/7_lean_RF_model/Figure_6A.pdf", plot = auc.box, width = 5, height = 5)
ggsave(filename = "output/7_lean_RF_model/Figure_6A.png", plot = auc.box, width = 5, height = 5, dpi = 400)


# 11. Classification when samples from a single cow are completely removed from the training dataset
####################################################################################################

# Register cluster 
cl <- makeCluster(10)
registerDoParallel(cl)

# "i" sample ID 
all.pred <- foreach(i = unique(samp.meta$CowN), .packages = "randomForest") %dopar%  {
    
    # Subset all samples from an animal
    c.id <- rownames(samp.meta)[samp.meta$CowN %in% i] 
    
    # Use a custom function to create training and validation datasets 
    all.d <- tr_and_valid_data(samp_valid = c.id,
                                 tr_subset_prop = 1,
                                 rf_data_in = rf.data.css.sig)
    
    # Create an RF model
    lean.rf <- randomForest(Shedder ~ ., 
             mtry = 28,
             data = all.d$training, 
             importance = TRUE,
             proximity = TRUE, 
             ntree = 15001) 
    
    # Predict status in training dataset 
       predict(lean.rf, all.d$validation[, !colnames(all.d$validation) %in% "Shedder"])
    
}

# 12. Plot classification results 
#################################

# Prepare data for plotting
# Unlist from object 
all.pred <- unlist(all.pred) 

# Combine vectors containing Cow Number, Cow ID, and predicted shedding from the unlisted object
c.pred.d <- cbind(sub("_.*", "", names(all.pred)), names(all.pred), as.character(all.pred))

# Adjust column names 
colnames(c.pred.d) <- c("CowN", "CowID", "Predicted_shedding")

# Convert WeightedScoreII column into categorical Shedding column 
samp.meta$Shedder <- ifelse(samp.meta$WeightedScoreII < 0.51, "Low", "High")

# Add CowID Column 
samp.meta$CowID <- rownames(samp.meta)

# Combine data about prediction and metadata 
pred.all.j <- left_join(as.data.frame(c.pred.d), samp.meta, by = "CowID")

# Add count column for plotting purposes 
pred.all.j$Count <- 1

# Evaluate if prediction is true and add results as a column
pred.all.j$Prediction <- pred.all.j$Predicted_shedding == pred.all.j$Shedder

# Arrange factors order (age in month decreasing)
pred.all.j$AgeMonth <- factor(pred.all.j$AgeMonth, levels = unique(pred.all.j$AgeMonth[order(as.numeric(pred.all.j$AgeMonth))]) )

# Plot and save graphs 
classific.cow.a <- ggplot(data = pred.all.j, aes(y = Count, x = CowN.x, fill = Prediction)) +
                     geom_bar(stat = "identity", color = "black") + 
                     facet_grid(~ Shedder, scales = "free") + 
                     theme_bw() + 
                     theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.title.x = element_blank()) +
                     theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())

classific.cow.b <- ggplot(data = pred.all.j, aes(y = Count, x = AgeMonth, fill = Prediction)) +
                         geom_bar(stat = "identity", color = "black") + 
                         facet_grid(CowN.x ~., switch = "y") + 
                         theme_bw() + 
                         theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
                         theme(strip.text.y = element_text(angle = 180)) + 
                         scale_y_continuous(breaks = NULL) 

ggsave(filename = "output/7_lean_RF_model/Figure_S7.png", plot = classific.cow.a, width = 10, height = 2.5, dpi = 400)
ggsave(filename = "output/7_lean_RF_model/Figure_S7.pdf", plot = classific.cow.a, width = 10, height = 2.5)

ggsave(filename = "output/7_lean_RF_model/Figure_S7.png", plot = classific.cow.b, width = 6, height = 5, dpi = 400)
ggsave(filename = "output/7_lean_RF_model/Figure_S7.pdf", plot = classific.cow.b, width = 6, height = 5)

# 13. Summarize information about prediction accuracy 
##################################################

# Create an empty list 
ind.cow.classific <- list()

# Summary via the table function
classific.sum <- data.frame(table(pred.all.j[,c(1, 19)]))

# Convert to wide format 
classific.sum <- spread(classific.sum, Prediction, Freq)

# Sum of prediction (n samples per cow)
classific.sum$N <- (as.numeric(classific.sum[,2]) + as.numeric(classific.sum[,3]))

# Percentage of false predictions 
classific.sum$False <- round(((as.numeric(classific.sum[, 2]) / classific.sum$N) * 100), 1) 

# Percentage of true predictions
classific.sum$True <- round(((as.numeric(classific.sum[, 3]) / classific.sum$N) * 100), 1) 

# Add shedding status 
classific.sum$Shedder <- pred.all.j$Shedder[! duplicated(pred.all.j$CowN.y)]

# Add summary table to summary list  
ind.cow.classific[["Summary_tabel"]] <- classific.sum

# Add summary about true prediction in the summary list 
ind.cow.classific[["Summary_true_prediciton"]] <- summary(classific.sum$True)

# Add information about CI to summary list 
ind.cow.classific[["ci_true_prediciton"]] <- ci(classific.sum$True)

