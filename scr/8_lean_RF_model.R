
#################
# "Lean" RF model
#################

# Set enviornment
#################

set.seed(49562)

source("scr/functions/general/load_abs_install_pkg.R")

load_abs_install_pkg(c("phyloseq", "tidyverse", "foreach", "doParallel", "randomForest"))

dir.create("output/8_lean_RF_model")


# Extract and format data for RF regression
# from phyloseq using a custom function 
#########################################

ps.tf.css.01 <- readRDS("output/3_filtering_mormalization/ps_tf2_css.RDS")

source("scr/functions/data_for_rf.R")

rf.data.css.01 <- data_for_rf(phyloseq = ps.tf.css.01, 
                              class.column = 'WeightedScoreII', 
                              remove.taxa.prev.less.than = 1, 
                              return.df = TRUE)

rf.data.css.01$WeightedScoreII <- as.numeric(as.character(rf.data.css.01$WeightedScoreII))


# Subset only significantly contributing taxa 
#############################################

sig.tax <- read.csv("output/7_discriminatory_tax/sig_contr_taxa.csv")

rf.data.css.sig <- rf.data.css.01[, as.character(sig.tax$X)]

rf.data.css.sig$WeightedScoreII <- rf.data.css.01$WeightedScoreII


# Create a list of datasets with different configuration of the shedding groups 
###############################################################################

score <- unique(rf.data.css.sig$WeightedScoreII)

score <- score[order(score, decreasing = TRUE)]

rfd.cat.l <- list()

for (s in 4:(length(score)-4)) {
    
    rfd.cat <- rf.data.css.sig
    
    colnames(rfd.cat)[colnames(rfd.cat) %in% "WeightedScoreII"] <- "Shedder"
    
    rfd.cat$Shedder <- as.factor(ifelse(rfd.cat$Shedder < score[s], "Low", "High"))
    
    rfd.cat.l[[s]] <- rfd.cat
    
}

# Make RF models 
################

cl <- makeCluster(16)
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

rf.plot.d <- c()

for (i in 1:length(rf.cat.l.res)) {
    
    rf.ind <- rf.cat.l.res[[i]]
    
    rf.plot.d <- rbind(rf.plot.d, cbind(rf.ind$confusion, rep(i, 2)))
}

rf.plot.d <- data.frame(rf.plot.d)

rf.plot.d$Class <- sub("\\.", "", rownames(rf.plot.d))

rf.plot.d$Class <- gsub("[0-9]", "", rf.plot.d$Class)

rf.plot.d$Split <- as.character(rep(score[4:(length(score)-4)], each = 2))


# Best category split plot 
##########################

lean.rf <- ggplot(rf.plot.d, aes(x=Split, y=class.error, color=Class)) + 

                                         geom_point(size=3) + 

                                         theme_bw() + 

                                         theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + 
                                        
                                         scale_color_manual(values = c("steelblue", "gold3"))

ggsave(filename = "output/8_lean_RF_model/best_cat_leanRF.png", plot = lean.rf, width = 5, height = 3.5, dpi = 300)

ggsave(filename = "output/8_lean_RF_model/best_cat_leanRF.pdf", plot = lean.rf, width = 5, height = 3.5)

# Determine best mtry and ntree for general rf model
####################################################

colnames(rf.data.css.sig)[colnames(rf.data.css.sig) %in% "WeightedScoreII"] <- "Shedder" 

rf.data.css.sig$Shedder <- as.factor(ifelse(as.numeric(as.character(rf.data.css.sig$Shedder)) < 0.51, "Low", "High"))

source("scr/functions/Tree_Mtry_Plot.R")

cl <- makeCluster(5)
registerDoParallel(cl)

tree.mtry.plot.all <- Tree_Mtry_Plot(data = rf.data.css.sig, 
                                     ntrees = c(7501, 10001, 15001), 
                                     start_val = round(sqrt(ncol(rf.data.css.sig)), 0), 
                                     stepF = 0.5, 
                                     class_colum = "Shedder", 
                                     ntimes = 5)

stopCluster(cl)


ggsave(filename = "output/8_lean_RF_model/lean_rf_mtry.png", plot = tree.mtry.plot.all, dpi = 400) 

ggsave(filename = "output/8_lean_RF_model/lean_rf_mtry.pdf", plot = tree.mtry.plot.all) 
