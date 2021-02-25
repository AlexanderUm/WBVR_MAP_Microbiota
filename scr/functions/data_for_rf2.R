##########################################################
#The funciton will extract otu table from the dataset,    #
#remove taxa with prevalence less than idetified,         #
#and add collumn for classification                       #
###########################################################
#* Taxa should be columns 

data_for_rf2 <- function(phyloseq, class.column, remove.taxa.prev.less.than) {
    
    #Extract otu table 
    otu.tab <- phyloseq@otu_table@.Data
    
    #Make prevalence otu table 
    otu.tab.prev <- otu.tab
    otu.tab.prev[otu.tab.prev > 0] <- 1 
    
    #Subset data
    otu.tab.f <- otu.tab[,colSums(otu.tab.prev) > remove.taxa.prev.less.than]
    
    #Add column 
     samp.data <- sample_data(phyloseq)
  
     otu.tab.comb <- cbind(otu.tab.f, samp.data[, class.column])

    #Return final object 
    return(otu.tab.comb)
} 