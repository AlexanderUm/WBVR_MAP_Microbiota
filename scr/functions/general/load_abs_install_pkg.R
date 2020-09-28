# The following function (load_abs_install_pkg) will load specified packages and install any missing form the list.
# BiocManager is used for installation

load_abs_install_pkg <- function(list.of.packages){
    
    if (!requireNamespace("BiocManager", quietly = TRUE))    
        install.packages("BiocManager")
    
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

    if(length(new.packages)) BiocManager::install(new.packages, suppressUpdates=TRUE)    
    
    lapply(list.of.packages, require, character.only = TRUE)
}
