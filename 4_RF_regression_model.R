##############################################################
# Prediction of shedding scores using random forest regression 
##############################################################


# Set seed 
set.seed(4957936)

# Load required libraries 
source("scr/functions/general/load_abs_install_pkg.R")

load_abs_install_pkg(c("phyloseq", "tidyverse", "foreach", "doParallel", "randomForest"))

# Create output directory 
dir.create("output/4_RF_regression_model")

# 1. Extract and format data for the RF regression
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



# 2. Find the optimal number of trees necessary for the model
#############################################################

# Make an RF regression model 
RF.tree <- randomForest(WeightedScoreII ~ . , 
                           data = rf.data.css.01, 
                           ntree = 15001)

# Visualize and save error rate in response to the number of used trees 
pdf("output/4_RF_regression_model/Figure_S4A.pdf")
plot(RF.tree)
dev.off()


# 3. Mtry parameter tuning 
###########################

# Load custom function for Mtry tuning and visualization 
source("scr/functions/Tree_Mtry_Plot.R")

# Make a computing cluster (number of cores)
cl <- makeCluster(10)

# Register the cluster 
registerDoParallel(cl)

# Visualize Mtry tunning attempts 
tree.mtry.plot.all <- Tree_Mtry_Plot(data = rf.data.css.01, 
                                     ntrees = 7501,
                                     start_val = ncol(rf.data.css.01/3), 
                                     stepF = 0.5, 
                                     class_colum = "WeightedScoreII", 
                                     ntimes = 3)

# Stop the computing cluster 
stopCluster(cl)

# Save the plot 
ggsave(plot = tree.mtry.plot.all, filename = "output/4_RF_regression_model/Figure_S4B.pdf")


# 4. Prediction of Shedding values could be regression mode of RF
#################################################################

# Load the function for creating matrices with randomly drawn not overlapping samples 
source("scr/functions/rand_draw_mat.R")

# Load the custom function for building RF models without a set of samples 
#           and consequently predict their values
source("scr/functions/rf_and_test.R")

# Make a computing cluster (number of cores)
cl <- makeCluster(10)

# Register the cluster 
registerDoParallel(cl)

# Create an empty list 
rf.reg.res <- list()

# Draw consequntly samples randomly substructed from from the datased
#      and predict their values. Samples will be drawn in five separate runs. 
for (i in 1:5) {

    # Create the matrix with randomly drawn not overlapping samples (10 samples per draw)
    rand.samp.reg <- rand_draw_mat(Samples_list = rownames(rf.data.css.01), Number_of_samp = 10)
    
    # Build an RF model without the randomly drawn set of samples and consequently predict their values. 
    res.for <- foreach(i=1:nrow(rand.samp.reg),.packages = "randomForest") %dopar% {
                rf_and_test(rf_data = rf.data.css.01, 
                samples_to_test = rand.samp.reg[i,], 
                n_samples_training = 230, 
                mtry = ncol(rf.data.css.01/3), 
                ntree = 7501, 
                variable_column = "WeightedScoreII", 
                regression_TorF = TRUE)}

    # Combine results 
    rf.reg.res[[i]] <- unlist(res.for) }

# Stop the computing cluster 
stopCluster(cl)

# Save the object 
save(rf.reg.res, file = "output/4_RF_regression_model/samples_reg_out.Rdata")


# 5. Extract and optimize data for visualization 
################################################

# Extract a part of metadata
p.reg.meta <- data.frame(ps.tf.css.01@sam_data[, c("CowN", "WeightedScoreII", "AgeMonth")])

# Adjust row names 
p.reg.meta$ID <- rownames(p.reg.meta)

# Create an empty dataframe for formatted data 
reg.pd <- data.frame()

# Extract data from the list and add metadata (loop)
# i is a list object number 
for(i in 1:length(rf.reg.res)) {
    
    # Extract individual object 
    reg.r <- data.frame(rf.reg.res[[i]])
    
    # Add an ID column containing rownames 
    reg.r$ID <- rownames(reg.r)
    
    # Add metadata and bind into a long dataframe 
    reg.pd <- rbind(reg.pd, left_join(reg.r, p.reg.meta, by="ID"))
    }

# Prepare data for the barchart (first layer) by 
#         leaving one row per cow 
reg.pd.bar <- reg.pd[!duplicated(reg.pd$CowN), ]

# Create vector with ordered Cow ID
ord.id <- reg.pd.bar$CowN[order(as.numeric(sub(",", ".", as.character(reg.pd.bar$WeightedScoreII))))]

# Arrange levels order for Cow ID 
reg.pd$CowN <- factor(reg.pd$CowN, levels = ord.id)

# Find mean value of predicted shedding scores per animal 
reg.pd.mean <- data.frame(aggregate(reg.pd$rf.reg.res..i.., list(reg.pd$CowN), mean))

# Adjust columns names 
colnames(reg.pd.mean) <- c("CowN", "PredV")


# 6. Visualize results of RF regression model prediction
########################################################

# Plot the results 
pred.val.plot <- ggplot() + 
    geom_bar(data = reg.pd[!duplicated(reg.pd$CowN), ], 
             aes(x = CowN, y = WeightedScoreII), 
             fill ="grey60", stat = "identity") + 
    geom_jitter(data = reg.pd, 
                aes(x = CowN, y = rf.reg.res..i.., color=AgeMonth), 
                width = 0.2) + 
    geom_smooth(data = reg.pd, 
                aes(x = as.numeric(CowN), y = rf.reg.res..i..), 
                method = "loess") + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) + 
    ylab("Weighted Score") + 
    xlab("Cow ID")

# Save plots
ggsave(filename = "output/4_RF_regression_model/Figure_S5.pdf", pred.val.plot, width = 7, height = 4)
ggsave(filename = "output/4_RF_regression_model/Figure_S5.png", pred.val.plot, width = 7, height = 4, dpi=400)


# 7. Summary of the RF regression results 
#########################################

# Summarize differences between predicted and actual value 
predv.dif.summ <- summary(abs(reg.pd$rf.reg.res..i.. - reg.pd$WeightedScoreII))

# Format the summary table 
predv.dif.summ <- as.data.frame(cbind(names(predv.dif.summ), 
                                      as.vector(predv.dif.summ)))

# Write summary table into file 
write.csv(predv.dif.summ, "output/4_RF_regression_model/predv_dif_summ.csv")



# 8. Correlation analysis between predicted and actual shedding values
######################################################################

# Prepare data for correlation analysis 
# Combine mean predicted values dataframe and data frame contained actual values 
cor.d <- inner_join(reg.pd.mean, reg.pd, by = "CowN")

# Remove rows with not unique cows names 
cor.d1 <- cor.d[!duplicated(cor.d$CowN), ]

# Correlation test between predicted and actual values.  
c.pred.score <- cor.test(cor.d1$PredV, cor.d1$WeightedScoreII)

# Correlation between differences in predicted and actual values and animals age.
c.sdiv.age <- cor.test(reg.pd$rf.reg.res..i.. - reg.pd$WeightedScoreII, reg.pd$AgeMonth)

# Correlation between differences in predicted and actual values and individual cows.
c.sdiv.cow <- cor.test(reg.pd$rf.reg.res..i.. - reg.pd$WeightedScoreII, as.numeric(reg.pd$CowN))

# Combine results of correlation analysis into a single dataframe.
all.cor <- rbind(c(c.pred.score$conf.int, c.pred.score$estimate, c.pred.score$p.value),
                    c(c.sdiv.age$conf.int, c.sdiv.age$estimate, c.sdiv.age$p.value),
                    c(c.sdiv.cow$conf.int, c.sdiv.cow$estimate, c.sdiv.cow$p.value))

# Adjust columns names in the combined table. 
colnames(all.cor) <- c("ci_low", "ci_high", "r(cor)", "p_val")

rownames(all.cor) <- c("Predicte & Real", "Age & Difference in values", "Cow ID & Difference in values")

# Write the combined table into a file. 
write.table(as.data.frame(all.cor), "output/4_RF_regression_model/correlations.txt")
