
##############################
# Investigate betta diversity 
##############################

# Load required libraries 
source("scr/functions/general/load_abs_install_pkg.R")

load_abs_install_pkg(c("phyloseq", "tidyverse", "vegan", "lubridate"))

dir.create("output/4_batta_diversity")


# Read in the phyloseq object 
ps.tf.css.01 <- readRDS("output/3_filtering_mormalization/ps_tf2_css.RDS")


###########################################
# Plot PCoA and test variables significance
###########################################

# Load function for data extraction 
source("scr/functions/pcoa_plot_data.R")


# Ordinate complite dataset 
###########################

# Extract data
all.samp <- pcoa_plot_data(phyloseq = ps.tf.css.01, 
                           factors = c("Year", "Diet"), 
                           vectors = c("AgeMonth", "WeightedScoreII"),
                           arrow_scale_cof = 3, 
                           ncore = 36)

# Plot ordination
all.samp.plot <- ggplot()+

     geom_point(data = all.samp$plot_df, aes(x = x, y = y, color = AgeMonth), size = 2) +

     theme_bw() +

     stat_ellipse(data = all.samp$plot_df, 
                  aes(x = x, y = y, fill = Diet), 
                  alpha = .1, size = 0.5, geom = "polygon") +

     geom_segment(data = all.samp$plot_vec, 
                  aes(x = 0, y = 0, xend = Dim1, yend = Dim2, lty = Vectors), 
                  arrow=arrow(length=unit(.2, "cm")))  + 

     xlab(label = all.samp$plot_var[1]) + 
     ylab(label = all.samp$plot_var[2])

# Save plots 
ggsave(filename = "output/4_batta_diversity/MDS_bray_all.pdf", plot = all.samp.plot, 
       device = "pdf", height = 7, width = 9)

ggsave(filename = "output/4_batta_diversity/MDS_bray_all.png", plot = all.samp.plot, 
       device = "png", height = 7, width = 9, dpi = 400)


# Ordination - solid diet 
#########################

ps.solid <- prune_samples(samples = ps.tf.css.01@sam_data$Diet %in% "Solid",
                           x = ps.tf.css.01)

# Extract data
solid.samp <- pcoa_plot_data(phyloseq = ps.solid, 
                           factors = c("Year"), 
                           vectors = c("AgeMonth", "WeightedScoreII"), 
                           arrow_scale_cof = 3, 
                           ncore = 36)


# Build and save the plot 
solid.samp.plot <- ggplot() +

     geom_point(data = solid.samp$plot_df,
                aes(x = x, y = y, color = AgeMonth), size=2) + 

     theme_bw() +

     geom_segment(data = solid.samp$plot_vec, 
                  aes(x = 0, y = 0, xend = Dim1, yend = Dim2, lty = Vectors), 
                  arrow=arrow(length=unit(.2, "cm"))) + 

     stat_ellipse(data = solid.samp$plot_df,
                  aes(x = x, y = y, fill=Year),
                  alpha=.1, size =0.5, geom="polygon") +

     xlab(label = solid.samp$plot_var[1]) +
     ylab(label = solid.samp$plot_var[2])  
     

ggsave(filename = "output/4_batta_diversity/MDS_bray_solid.pdf", plot = solid.samp.plot, 
       device = "pdf", height = 5, width = 8)

ggsave(filename = "output/4_batta_diversity/MDS_bray_solid.png", plot = solid.samp.plot, 
       device = "png", height = 5, width = 8, dpi = 400)
