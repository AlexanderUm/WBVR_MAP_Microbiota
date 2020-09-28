###########################################################
#The funciton will extract otu table from the dataset,    #
#remove taxa with prevalence less than idetified,         #
#and add collumn for classification                       #
###########################################################
#* Taxa should be columns 

data_for_rf <- function(phyloseq, class.column, remove.taxa.prev.less.than, return.df) {
    
    #Extract otu table 
    otu.tab <- phyloseq@otu_table@.Data
    
    #Make prevalence otu table 
    otu.tab.prev <- otu.tab
    otu.tab.prev[otu.tab.prev > 0] <- 1 
    
    #Subset data
    otu.tab.f <- otu.tab[,colSums(otu.tab.prev) > remove.taxa.prev.less.than]
    
    #Add column 
     samp.data <- data.frame(phyloseq@sam_data@.Data)
     colnames(samp.data) <- phyloseq@sam_data@names
     otu.tab.comb <- cbind(otu.tab.f, as.character(samp.data[, class.column]))
     colnames(otu.tab.comb) <- c(colnames(otu.tab.f), class.column)
    
    if (return.df == TRUE) {
        otu.tab.comb <- data.frame(otu.tab.comb)
        for (i in colnames(otu.tab.f)) {
        otu.tab.comb[,i] <- as.numeric(as.character(otu.tab.comb[,i]))}
    }
    
    #Return final object 
    return(otu.tab.comb)
} 
