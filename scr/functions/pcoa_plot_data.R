# Following function will ordinate supplied data, fit variables, 
# test significance of variables, and extract data for plotting 
# * bray dissimularity matrix is used for ordination 
################################################################

pcoa_plot_data <- function(phyloseq, factors, vectors, arrow_scale_cof, ncore, nperm) {
    
    # Extract data
    ###############
    
    # Extract OTU table 
    asv.tab0 <- otu_table(phyloseq)
    
     # Extract sample table 
    meta.d <- data.frame(sample_data(phyloseq))
    
    # Select variables to test 
    
    meta.d.fit <- mutate_all(meta.d[, c(factors, vectors)], as.character)
    
    meta.d.fit <- mutate_at(meta.d.fit, factors, as.factor)
    
    meta.d.fit <- mutate_at(meta.d.fit, vectors, as.numeric)
    
    
    # Construct ordination
    #######################
    
    # Create bray dissimularity matrix 
    asv.tab0.bray <- vegdist(asv.tab0, method = "bray")
    
    # MDS on dissimularity matrix = PCoA
    asv.tab0.bray.mds <- cmdscale(asv.tab0.bray, eig = TRUE, x.ret = TRUE)
    
     # Fit variables into ordination
    fit <- envfit(ord = asv.tab0.bray.mds, env = meta.d.fit)
    
    # Prepare data for plotting
    ###########################
    
    # Extract % of variations captured by axes 
    axis.cover <- round(asv.tab0.bray.mds$eig/sum(asv.tab0.bray.mds$eig)*100, 1)
    
    axis.cover[1] <- paste0("[ ", axis.cover[1], "%",  " ]" )
    
    axis.cover[2] <- paste0("[ ", axis.cover[2], "%",  " ]" )
    
    # Extract vector coordinates 
    arrow <- data.frame(fit$vectors$arrows/arrow_scale_cof, 
                        R = fit$vectors$r, 
                        P = fit$vectors$pvals)

    arrow$Vectors <- rownames(arrow)
    
    # Extract centroids coordinates 
    text.df <- data.frame(fit$factors$centroids)
    
    # Collect data for the plot 
    asv.mds.val <- asv.tab0.bray.mds$points

    asv.bray.mds.df <- data.frame(Samples = rownames(asv.mds.val), 
                              x = asv.mds.val[,1], 
                              y = asv.mds.val[,2],
                              meta.d.fit)   
    
    # Test significance of varialbes fitting 
    ########################################
    
    # Define variables if missing 
    if (missing(ncore)) { ncore = 1 }
    
    if (missing(nperm)) { nperm = 999}
    
    # Test significance 
    adonis.out <- adonis2(asv.tab0 ~ ., 
                           data = meta.d.fit,
                           permutations = nperm, 
                           parallel = ncore,
                           by = "terms")
    
    # Collect objects in a list 
    ret.obj <- list("adonis2_test" = adonis.out, 
                    "plot_df" = asv.bray.mds.df, 
                    "plot_vec" = arrow, 
                    "plot_var" = axis.cover)
    
    return(ret.obj)

}