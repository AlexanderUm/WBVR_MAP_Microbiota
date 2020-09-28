
setwd('/home/umane001/WBVR_MAP_Microbiota')

#Load required libraries 
#########################

source("scr/functions/general/load_abs_install_pkg.R")

load_abs_install_pkg(c("phyloseq", "tidyverse", "foreach", "doParallel", "randomForest"))

dir.create("output/9_time_influence")

###########################
# Prepare data for RF model 
###########################

# Creating an empty list where all datasets will be stored
data.sets.list <- list()

# Creating an empty list where corresponding metadata will be stored
meta.sets.list <- list()

# Prepare metadata 
##################

ps.tf.css.01 <- readRDS("output/objects/phyloseq/ps_tf_css_011.RDS")

ps.tf.meta <- data.frame(sample_data(ps.tf.css.01))

ps.tf.meta$YearMonth <- paste(ps.tf.meta$Year_shed, ps.tf.meta$Month_shed, sep="_")


# RF data full model / all samples 
###################################

source("scr/functions/data_for_rf.R")


rf.data.css <- data_for_rf(phyloseq = ps.tf.css.01, 
                              class.column = 'WeightedScoreII', 
                              remove.taxa.prev.less.than = 1, 
                              return.df = TRUE)

rf.data.css$WeightedScoreII <- as.numeric(as.character(rf.data.css$WeightedScoreII))

colnames(rf.data.css)[colnames(rf.data.css) %in% "WeightedScoreII"] <- "Shedder"

rf.data.css$Shedder <- as.factor(ifelse(rf.data.css$Shedder < 0.51, "Low", "High"))

data.sets.list[["atax.asamp"]] <- rf.data.css
meta.sets.list[["m.atax.asamp"]] <- ps.tf.meta

# RF data lean model / all samples 
###################################

sig.tax <- read.csv("output/plots/7_discriminatory_tax/sig_contr_taxa.csv")

data.sig <- rf.data.css[, as.character(sig.tax$X)]

data.sig$Shedder <- rf.data.css$Shedder

data.sets.list[["ltax.asamp"]] <- data.sig
meta.sets.list[["m.ltax.asamp"]] <- ps.tf.meta

# RF data full model / long lived cows 
######################################

meta.longl <- ps.tf.meta[ps.tf.meta$CowN_shed %in% c("1349", "1350", "1351", "1355", "1356", "1357", "1359", 
                                       "1360", "1362", "1363", "1364", "1367"), ] 

rf.data.longl <- rf.data.css[rownames(meta.longl), ]

shedd.longl <- rf.data.longl$Shedder

rf.data.longl <- rf.data.longl[,!colnames(rf.data.longl) %in% "Shedder"]

rf.data.longl <- rf.data.longl[,!colSums(rf.data.longl) == 0]

rf.data.longl$Shedder <- shedd.longl

data.sets.list[["atax.lsamp"]] <- rf.data.longl
meta.sets.list[["m.atax.lsamp"]] <- meta.longl

# RF data lean model / long lived cows 
######################################

data.sig.longl <- rf.data.longl[, colnames(rf.data.longl) %in% as.character(sig.tax$X)]

data.sig.longl <- data.sig.longl[,!colSums(data.sig.longl) == 0]

data.sig.longl$Shedder <- rf.data.longl$Shedder

data.sets.list[["ltax.lsamp"]] <- data.sig.longl
meta.sets.list[["m.ltax.lsamp"]] <- meta.longl

########################################
# Run RF analysis on time based datasets
########################################

#Load functions 
###############
source("scr/functions/rf_by_time.R")

source("scr/functions/rf_by_size.R")


save_plot <- function(directory, plot, w, h, pnam){
    
    dirnam <- paste0(directory, "/", gsub("\\.", "_", pnam))
    
    ggsave(filename = paste0(dirnam, ".png"), plot = plot, width = w, height = h, dpi = 400)
    
    ggsave(filename = paste0(dirnam, ".pdf"), plot = plot, width = w, height = h)
    
    plot.c <- plot + coord_cartesian(xlim = c(0, 1), ylim=c(0, 1), expand = TRUE)
    
    ggsave(filename = paste0(dirnam, "_exp", ".png"), plot = plot.c, width = w, height = h, dpi = 400)
    
    ggsave(filename = paste0(dirnam, "_exp", ".pdf"), plot = plot.c, width = w, height = h)
    
    write.csv(plot[["data"]], paste0(dirnam, ".csv"))
}


# Prepare support objects 
mtr.v <- c(285, 9, 285, 9)

n.tree <- c(151, 101, 151, 101)

cat.colors <- c("steelblue", "gold")

memory.obj <- ls()

###################################
# RF time-slice model 
# All taxa / all samples 
###################################
 
for (d in 1:length(data.sets.list)) {

    cl <- makeCluster(56)
    registerDoParallel(cl)  
    
    rand <- rf_by_size(data = data.sets.list[[d]], 
                sub_percent = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 ), 
                class_colum = "Shedder", 
                class_colors = cat.colors, 
                mtr = mtr.v[d], 
                ntree = n.tree[d], 
                nrep = 2)
    
    new.o <- paste0("rand.", names(data.sets.list)[d])
    
    assign(new.o, rand)
    
    save_plot(directory = "output/9_time_influence/plots", plot = get(new.o), w = 10, h = 5.5, pnam = new.o)

stopCluster(cl)
obj.rem <- ls()[!ls() %in% c(memory.obj, "d")]
rm(list=obj.rem)
gc()

    
cl <- makeCluster(56)
registerDoParallel(cl)
    
from.early <- rf_by_time(data_for_rf = data.sets.list[[d]],
           metadata_for_rf = meta.sets.list[[d]],
           class_column = "Shedder",
           class_colors = cat.colors,
           slice_column = "YearMonth",
           slice_steps = c('1999_02', '1999_05', '1999_09', '2000_04', 
                           '2000_08', '2001_01', '2001_05', '2001_08', '2001_09', '2001_12', 
                          '2002_03', '2002_05', '2002_09', '2002_10', '2003_02'), 
           mtr = mtr.v[d],
           ntree = n.tree[d], 
           nrep = 2)
    
    new.o <- paste0("from.early.", names(data.sets.list)[d])
    
    assign(new.o, from.early)
    
    save_plot(directory = "output/9_time_influence/plots", plot = get(new.o), w = 10, h = 5.5, pnam = new.o)
    
stopCluster(cl)
obj.rem <- ls()[!ls() %in% c(memory.obj, "d")]
rm(list=obj.rem)
gc()
    
    
    
cl <- makeCluster(56)
registerDoParallel(cl)
    
from.late <- rf_by_time(data_for_rf = data.sets.list[[d]],
           metadata_for_rf = meta.sets.list[[d]],
           class_column = "Shedder",
           class_colors = cat.colors,
           slice_column = "YearMonth",
           slice_steps = c('2003_05', '2003_03', '2003_02', '2002_10', '2002_09', '2002_05', 
                           '2002_03', '2001_12', '2001_09', '2001_08', 
                          '2001_05', '2001_01', '2000_08', '2000_04', '1999_09'),  
           mtr = mtr.v[d],
           ntree = n.tree[d], 
           nrep = 2)

    new.o <- paste0("from.late.", names(data.sets.list)[d])
    
    assign(new.o, from.late)
    
    save_plot(directory = "output/9_time_influence/plots", plot = get(new.o), w = 10, h = 5.5, pnam = new.o)
 
stopCluster(cl)
obj.rem <- ls()[!ls() %in% c(memory.obj, "d")]
rm(list=obj.rem)
gc()

}
