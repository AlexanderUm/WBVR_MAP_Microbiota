
#Function: Make matrix with n randomly drown from data set samples per row
#          Samples are not repeating in each row
##########################################################################

rand_draw_mat <- function(Samples_list, Number_of_samp) {
    
    Samples_list_red <- Samples_list
    
    s.out <- c()
    
    s.out.all <- c()
    
    for (i in 1:round(length(Samples_list)/Number_of_samp, 0))  {
        
        Samples_list_red <- Samples_list_red[!Samples_list_red %in% s.out]
        
        if (length(Samples_list_red) > length(s.out)) {
            
             s.out <- sample(Samples_list_red, Number_of_samp)
    
             s.out.all <- rbind(s.out.all, s.out)
             }
        
        }
    
    return(s.out.all)
}