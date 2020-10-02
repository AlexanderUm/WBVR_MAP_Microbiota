
##############################################################
# 6. General random  forest model 
# - Find the best assign shedding class based on accuracy 
#   of classification. 
##############################################################

# Load required libraries 
#########################

set.seed(4597234)

source("scr/functions/general/load_abs_install_pkg.R")

load_abs_install_pkg(c("phyloseq", "tidyverse", "foreach", "doParallel",
                       "randomForest", "rfPermute", "ComplexHeatmap"))

dir.create("output/6_RF_general_model")

# Extract and format data for RF model
######################################

ps.tf.css.01 <- readRDS("output/3_filtering_mormalization/ps_tf2_css.RDS")

source("scr/functions/data_for_rf.R")

rf.data.css.01 <- data_for_rf(phyloseq = ps.tf.css.01, 
                              class.column = 'WeightedScoreII', 
                              remove.taxa.prev.less.than = 1, 
                              return.df = TRUE)

rf.data.css.01$WeightedScoreII <- as.numeric(as.character(rf.data.css.01$WeightedScoreII))


# Create a list of datasets with different configuration of shedding groups 
###########################################################################

# Extract and order shedding score
score <- unique(rf.data.css.01$WeightedScoreII)

score <- score[order(score, decreasing = TRUE)]

# Create datasets with different configuration of High and Low groups
rfd.cat.l <- list()

for (s in 4:(length(score)-4)) {
    
    rfd.cat <- rf.data.css.01
    
    colnames(rfd.cat)[colnames(rfd.cat) %in% "WeightedScoreII"] <- "Shedder"
    
    rfd.cat$Shedder <- as.factor(ifelse(rfd.cat$Shedder < score[s], "Low", "High"))
    
    rfd.cat.l[[s]] <- rfd.cat
}


# Make RF models 
################

cl <- makeCluster(24)

registerDoParallel(cl)

rf.cat.l.res <- foreach(p=4:length(rfd.cat.l),.packages = "randomForest") %dopar% {
    
                randomForest(Shedder ~ ., 
                             data=rfd.cat.l[[p]], 
                             importance=TRUE,
                             proximity=TRUE, 
                             ntree=7501) 
                }

stopCluster(cl)


# Extract and format data for plotting 
######################################

# Extract data from the list 
rf.plot.d <- c()

for (i in 1:length(rf.cat.l.res)) {
    
    rf.ind <- rf.cat.l.res[[i]]
    
    rf.plot.d <- rbind(rf.plot.d, cbind(rf.ind$confusion, rep(i, 2)))
}

# Format data for plotting 
rf.plot.d <- data.frame(rf.plot.d)

rf.plot.d$Class <- sub("\\.", "", rownames(rf.plot.d))

rf.plot.d$Class <- gsub("[0-9]", "", rf.plot.d$Class)

rf.plot.d$Split <- as.character(rep(score[4:(length(score)-4)], each = 2))


# Plot best category             
split.rf <- ggplot(rf.plot.d, aes(x=Split, y=class.error, color=Class)) + 

                geom_point(size=3) + 

                theme_bw() + 

                theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + 
                                            
                scale_color_manual(values = c("steelblue", "gold3"))

ggsave(filename = "output/6_RF_general_model/best_cat_rf.png", 
       plot = split.rf, width = 5, height = 3.5, dpi = 400)

ggsave(filename = "output/6_RF_general_model/best_cat_rf.pdf", 
       plot = split.rf, width = 5, height = 3.5)


#################################################
# Best possible RF model for the complite dataset
#################################################


# Prepare data for RF model 
###########################

ps.tf.css.01 <- readRDS("output/3_filtering_mormalization/ps_tf2_css.RDS")

source("scr/functions/data_for_rf.R")

rf.data.css <- data_for_rf(phyloseq = ps.tf.css.01, 
                              class.column = 'WeightedScoreII', 
                              remove.taxa.prev.less.than = 1, 
                              return.df = TRUE)

rf.data.css$WeightedScoreII <- as.numeric(as.character(rf.data.css$WeightedScoreII))

colnames(rf.data.css)[colnames(rf.data.css) %in% "WeightedScoreII"] <- "Shedder"

rf.data.css$Shedder <- as.factor(ifelse(rf.data.css$Shedder < 0.51, "Low", "High"))


# Determine best mtry and ntree for general rf model
####################################################

source("scr/functions/Tree_Mtry_Plot.R")

cl <- makeCluster(36)

registerDoParallel(cl)

tree.mtry.plot.all <- Tree_Mtry_Plot(data = rf.data.css, 
                                     ntrees = c(7001, 10001, 15001), 
                                     start_val = round(sqrt(ncol(rf.data.css)), 0), 
                                     stepF = 0.5, 
                                     class_colum = "Shedder", 
                                     ntimes = 5)

stopCluster(cl)


ggsave(filename = "output/6_RF_general_model/general_rf_mtry.png", plot = tree.mtry.plot.all, 
       width = 5, height = 4,  dpi = 400) 

ggsave(filename = "output/6_RF_general_model/general_rf_mtry.pdf", plot = tree.mtry.plot.all, 
       width = 5, height = 4) 

