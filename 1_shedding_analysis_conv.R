##################################################
# Analysis and visualisation of shedding patterns 
##################################################


# Set environment
#################

set.seed(3432)

source("scr/functions/general/load_abs_install_pkg.R")

load_abs_install_pkg(c("phyloseq", "tidyverse", "plyr"))

dir.create("output/1_shedding_analysis")


# 1. Read and format shedding data for analysis
###############################################

shed.data <- read.csv("data/shedding_data.csv", 
                      stringsAsFactors = FALSE, 
                      na.strings=c("","NA"))

# Format columns names 
colnames(shed.data) <- sub("X", "C", colnames(shed.data))

# Format date 
shed.data$SDATE <- gsub("-", "_", shed.data$SDATE)

# Melt data to the long format 
shed.data.l <- gather(shed.data, CowID, Score, C1348:C1367)

# Remove data points without shedding values
shed.data.l <- shed.data.l[shed.data.l$Score %in% c("0", "+", "++", "+++"), ]

# Save formatted and trimmed data
write.csv(shed.data.l, "output/1_shedding_analysis/shedding_data_long.csv")


# 2. Calculate shedding scores
##############################

# Convert shedding score to numeric values
shed.data.l$ScoreNum <- as.numeric(as.character(mapvalues(shed.data.l$Score, 
                                                          c("0", "+", "++", "+++"), 
                                                          c(0, 1, 2, 3))))

# Calculate cumulative shedding score per animal
samp.score <- aggregate(shed.data.l$ScoreNum, 
                        by=list(Category=shed.data.l$CowID), FUN=sum)

# Calculate weighted shedding scores
w.score <- round(samp.score[,"x"]/table(shed.data.l$CowID), 2)





# 3. Calculate weighted shedding score for the first 43 time point 
##################################################################

# Select first 43 time points 
retain.points <- shed.data$SDATE[1:43]

# Subset shedding data leaving only samples from the first 43 data points
shed.data.l2 <- shed.data.l[shed.data.l$SDATE %in% retain.points, ]

# Calculate sum of shedding scores per cow 
samp.score2 <- aggregate(shed.data.l2$ScoreNum, 
                        by = list(Category=shed.data.l2$CowID), FUN = sum)

# Calculate relative shedding score per cow (shedding scores sum / by number of samples)
w.score2 <- round(samp.score2[,"x"] / table(shed.data.l2$CowID), 2)


# 4. Combine information about shedding and write it into a file 
################################################################

# Bind  Cumulative shedding score, Weighted shedding score overall, and 
#       Weighted shedding score for the first 43 data points into a dataframe 
all.scores <- cbind(samp.score, as.vector(w.score), as.vector(w.score2))

# Add column names 
colnames(all.scores) <- c("CowN", "CumulativeScore", "WeightedScoreI", "WeightedScoreII")

# Write the dataframe into a file 
write.csv(all.scores, "output/1_shedding_analysis/shedding_scores.csv")



# 5. Add the information about shedding to general metadata
###########################################################

# Read in metadata file 
s.meta <- read.csv("data/metadata/samples_metadata_f.csv")

# Add shedding data to the metadata 
s.meta.comb <- left_join(s.meta, all.scores, by="CowN")

# Write updated metadata into the file 
write.csv(s.meta.comb[,-1], "output/1_shedding_analysis/metadata_f.csv")



# 6. Correlation between w.score and w.score2 
#############################################

# Select only cows that survive till the end of experiment  
long.lived.animals <- c("C1349", "C1350", "C1351", "C1355", "C1356", "C1357", 
                        "C1359", "C1360", "C1362", "C1363", "C1364", "C1367")

# Correlate Overall Weighted shedding scores and 
#           Weighted shedding scores for the first 43 time points 
#           using Pearson correlation 
corr.w.score <- cor.test(w.score2[long.lived.animals], 
                round(w.score[long.lived.animals], 2), 
                method = "pearson") 

# Select only P value and estimate from the correlation object 
corr.w.score.out <- c(corr.w.score$p.value, corr.w.score$estimate)

# Adjust names 
names(corr.w.score.out) <- c("Pval", "Estimate")

# Write the correlation analysis results into the file 
write.csv(corr.w.score.out, "output/1_shedding_analysis/corr_w_scores.csv")


# 7. Visualize Weighted shedding scores for the first 43 time points
####################################################################

# Convert data into a dataframe 
w.score2.dp <- data.frame(w.score2)

# Order cow IDs by shedding score 
w.score2.dp$Var1 <- factor(w.score2.dp$Var1, levels = w.score2.dp$Var1[order(w.score2.dp$Freq)])

# Create a column that reflects animals' life spend 
w.score2.dp$Life_spend <- ifelse(w.score2.dp$Var1 %in% long.lived.animals, "Normal", "Early Culled")

# Plot Weighted shedding scores for the first 43 time points
w.score2.p <- ggplot(w.score2.dp, aes(x = Var1, y = Freq, fill = Life_spend)) + 
  geom_bar(stat = "identity") + 
  theme_bw() +
  geom_hline(yintercept = mean(w.score2.dp$Freq), color="black") +
  geom_hline(yintercept = mean(w.score2.dp$Freq) + sd(w.score2.dp$Freq), color="black", linetype = "longdash") +
  geom_hline(yintercept = mean(w.score2.dp$Freq) - sd(w.score2.dp$Freq), color="black", linetype = "longdash") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  guides(fill = guide_legend(title="Animals Life Spend")) + 
  scale_fill_brewer(palette = "Accent") +
  ylab("Weighted score") + 
  xlab("Cow ID") + 
  coord_flip()

# Save plot into files as png and pdf formats 
ggsave(filename = "output/1_shedding_analysis/Figure_2A.pdf", w.score2.p, width = 4, height = 4)
ggsave(filename = "output/1_shedding_analysis/Figure_2A.png", w.score2.p, width = 4, height = 4, dpi = 400)


# 8. Visualize Overall Weighted shedding scores and 
#              Weighted shedding scores for the first 43 time points
#              together 
####################################################################

# Convert Overall Weighted shedding scores into a dataframe 
w.score1.dp <- data.frame(w.score)

# Create a column that reflects animals' life spend
w.score1.dp$Life_spend <- ifelse(w.score1.dp$Var1 %in% long.lived.animals, "Normal", "Early Culled")

# Add Score_type column to Overall Weighted shedding scores dataframe 
w.score1.dp$Score_type <- "Overall"

# Calculate and add shedding scores mean to Overall Weighted shedding scores dataframe 
w.score1.dp$Mean <- mean(w.score1.dp$Freq)

# Calculate and add columns with mean + and - standard deviation for
#          for Overall Weighted shedding scores dataframe
w.score1.dp$SD1 <- (mean(w.score1.dp$Freq) + sd(w.score1.dp$Freq))
w.score1.dp$SD2 <- (mean(w.score1.dp$Freq) - sd(w.score1.dp$Freq))

# Add Score_type column to Weighted shedding scores for the first 43 time points dataframe
w.score2.dp$Score_type <- "Normalized"

# Calculate and add shedding scores mean to Weighted shedding scores for the first 43 time points dataframe 
w.score2.dp$Mean <- mean(w.score2.dp$Freq)

# CCalculate and add columns with mean + and - standard deviation 
#          for Weighted shedding scores for the first 43 time points dataframe 
w.score2.dp$SD1 <- (mean(w.score2.dp$Freq) + sd(w.score2.dp$Freq))
w.score2.dp$SD2 <- (mean(w.score2.dp$Freq) - sd(w.score2.dp$Freq))

# Combine Overall Weighted shedding scores dataframe and 
#         Weighted shedding scores for the first 43 time points dataframe 
#         into a single long format dataframe 
w.score.comb.dp <- rbind(w.score2.dp, w.score1.dp)

# Plot Weighted shedding scores for the first 43 time points and 
#      Overall Weighted shedding scores together 
w.score.comb <- ggplot(w.score.comb.dp, aes(x=Var1, y=Freq, fill = Life_spend)) + 
  geom_bar(stat="identity") + 
  theme_bw() + 
  geom_hline(aes(yintercept = Mean), color="black") +
  geom_hline(aes(yintercept = SD1), color="black", linetype = "longdash") +
  geom_hline(aes(yintercept = SD2), color="black", linetype = "longdash") + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) + 
  guides(fill = guide_legend(title="Animals Life Spend")) + 
  scale_fill_brewer(palette = "Accent") +
  ylab("Weighted score") + 
  xlab("Cow ID") + 
  facet_grid(.~ Score_type) +
  coord_flip()

# Save plot into files as png and pdf formats 
ggsave(filename = "output/1_shedding_analysis/Figure_S1.pdf", w.score.comb, width = 5, height = 3.5)
ggsave(filename = "output/1_shedding_analysis/Figure_S1.png", w.score.comb, width = 5, height = 3.5, dpi = 400)


# 9. Calculate the increase in shedding score over time
########################################################

# Calculate age of cows in days using 1999-01-31 as a starting date 
tdif <- difftime(as.character(gsub("_", "-", shed.data.l$SDATE)), "1999-01-31", units = "days")

# Calculate age of cows in month by dividing age in days
#           by average number of days in a month
shed.data.l$AgeMonth <- round(as.numeric(tdif) / (365/12), 0)

# Calculate age of cows in quarters (3, 6, 9, 12) by taking the ceiling 
#           of age in month divided by 3 and then multiplying by 3 
shed.data.l$AgeMonthQ <- ceiling(shed.data.l$AgeMonth / 3) * 3

# Make an empty vector 
cshed.v <- c()

# Calculate cumulative shedding score per time poning (in a loop)
# i is cow ID 
for (i in unique(shed.data.l$CowID)) {
    
    # Subset shedding data for an individual cow 
    cowdat <- shed.data.l[shed.data.l$CowID %in% i, ]
    
    # Sum shedding scores for each quartal 
    sumQ <- tapply(cowdat$ScoreNum, cowdat$AgeMonthQ, FUN = sum)
    
    # Create an empty vector for cumulative scores
    sumQcumVec <- c()
    
    # Nested loop - cumulative score calculation 
    # n is quartal 
    for (n in 1:length(sumQ)) {
        
        # If it is first quartal create object sumQcum equal 
        #          first shedding score for the first quartal 
        if (n == 1) {sumQcum <- sumQ[1]}
        
        # Esle sumQcum equal previous sumQcum plus current 
        else {sumQcum <- sumQcum + sumQ[n]}
        
    # Combine cumulative scores in a vector      
    sumQcumVec <- c(sumQcumVec, sumQcum)
        
    }
    
    # Combine vector of cumulative scores,  
    #         vector of sums of shedding scores for each quartal, and 
    #         vector containing age in quartas into a dataframe 
    CumSum <- data.frame(cbind(sumQcumVec, sumQ, unique(cowdat$AgeMonthQ)))
    
    # Add cow ID 
    CumSum$CowN <- i 
    
    # Combine every entry per cow into a single long dataframe 
    cshed.v <- rbind(cshed.v, CumSum)
    
  } 


# 10. Visualisation of the increase in shedding score over time
####################################################

# Adjust column names for the Increase in shedding score over time dataframe 
colnames(cshed.v) <- c("Cumulative_SIS", "W_SIS", "Age_Month", "CowN")

# Adjust column names for the Weighted shedding scores for the first 43 time points dataframe 
colnames(w.score2.dp) <- c("CowN", "Freq", "Life_spend")

# Add the Weighted shedding scores for the first 43 time points dataframe to 
#     the Increase in shedding score over time dataframe 
cshed.v2 <- left_join(cshed.v, w.score2.dp[,1:3], by = "CowN")

# Set levels order for the Life_spend column 
cshed.v2$Life_spend <- factor(cshed.v2$Life_spend, levels = c("Normal", "Early Culled"))

# Create column shedding to indicate "High" and "Low" shedders 
#     [based on the future analysis cows that have 
#      Weighted shedding scores for the first 43 time points are "Low" and above are "High" shedders]
cshed.v2$Shedding <- ifelse(cshed.v2$Freq < 0.51, "Low", "High")

# Plot the increase in shedding score over time
cum.shed.p <- ggplot(cshed.v2, aes(x = Age_Month, y = Cumulative_SIS, group = CowN, color = Shedding)) + 
        geom_point(size = 1.2) + 
        geom_line( size = 0.5) + 
        facet_grid(Life_spend ~ .) +
        theme_bw() + 
        ylab("Cumulative shedding score") + 
        xlab("Age (Month)")  + 
        scale_color_manual(values = c("steelblue", "gold3")) + 
        scale_x_continuous(breaks = seq(0, 60, by = 6))

# Save plot into files as png and pdf formats 
ggsave(filename = "output/1_shedding_analysis/Figure_2B.pdf", cum.shed.p, width = 6, height = 4)
ggsave(filename = "output/1_shedding_analysis/Figure_2B.png", cum.shed.p, width = 6, height = 4, dpi = 400)

