
##################################################
# Analysis and visualisation of shedding patterns 
##################################################


# Set enviornment
#################

set.seed(3432)

source("scr/functions/general/load_abs_install_pkg.R")

load_abs_install_pkg(c("phyloseq", "tidyverse", "cowplot", "plyr"))

dir.create("output/1_shedding_analysis")

# Read and format shedding data for analysis
############################################

shed.data <- read.csv("data/shedding_data.csv", 
                      stringsAsFactors = FALSE, 
                      na.strings=c("","NA"))

colnames(shed.data) <- sub("X", "C", colnames(shed.data))

shed.data$SDATE <- gsub("-", "_", shed.data$SDATE)

# Melt data to the long format 
shed.data.l <- gather(shed.data, CowID, Score, C1348:C1367)

# Remove data points without shedding values
shed.data.l <- shed.data.l[shed.data.l$Score %in% c("0", "+", "++", "+++"), ]

# Save formated and trimed data
write.csv(shed.data.l, "output/1_shedding_analysis/shedding_data_long.csv")

# Calculate shedding scores
###########################

# Convert shedding score to numeric values
shed.data.l$ScoreNum <- as.numeric(as.character(mapvalues(shed.data.l$Score, 
                                                          c("0", "+", "++", "+++"), 
                                                          c(0, 1, 2, 3))))

# Calculate cumulative shedding score per animal
samp.score <- aggregate(shed.data.l$ScoreNum, 
                        by=list(Category=shed.data.l$CowID), FUN=sum)


# Calculate weighted shedding scores
w.score <- round(samp.score[,"x"]/table(shed.data.l$CowID), 2)



# Calculate weighted shedding score for the first 43 time point 
retain.points <- shed.data$SDATE[1:43]

shed.data.l2 <- shed.data.l[shed.data.l$SDATE %in% retain.points, ]

samp.score2 <- aggregate(shed.data.l2$ScoreNum, 
                        by=list(Category=shed.data.l2$CowID), FUN=sum)

w.score2 <- round(samp.score2[,"x"]/table(shed.data.l2$CowID), 2)


# Combine information about shedding and save it 
all.scores <- cbind(samp.score, as.vector(w.score), as.vector(w.score2))

colnames(all.scores) <- c("CowN", "CumulativeScore", "WeightedScoreI", "WeightedScoreII")

write.csv(all.scores, "output/1_shedding_analysis/shedding_scores.csv")


# Add information about shedding to general metadata
s.meta <- read.csv("data/metadata/samples_metadata_f.csv")

s.meta.comb <- left_join(s.meta, all.scores, by="CowN")

write.csv(s.meta.comb[,-1], "output/1_shedding_analysis/metadata_f.csv")


# Correlation w.score and w.score2 
###################################

# Select only cows that survive till the end of experiment  
long.lived.animals <- c("C1349", "C1350", "C1351", "C1355", "C1356", "C1357", 
                        "C1359", "C1360", "C1362", "C1363", "C1364", "C1367")

# Crrelation analysis
corr.w.score <- cor.test(w.score2[long.lived.animals], 
                round(w.score[long.lived.animals], 2), 
                method = "pearson") 

corr.w.score.out <- c(corr.w.score$p.value, corr.w.score$estimate)

names(corr.w.score.out) <- c("Pval", "Estimate")

write.csv(corr.w.score.out, "output/1_shedding_analysis/corr_w_scores.csv")


# Visualize shedding scores 
############################

# Plot w.score2
w.score2.dp <- data.frame(w.score2)

w.score2.dp$Var1 <- factor(w.score2.dp$Var1, levels=w.score2.dp$Var1[order(w.score2.dp$Freq)])

w.score2.p <- ggplot(w.score2.dp, aes(x=Var1, y=Freq)) + 
  geom_bar(stat="identity") + 
  theme_bw() +
  geom_hline(yintercept = mean(w.score2.dp$Freq), color="blue") +
  geom_hline(yintercept = mean(w.score2.dp$Freq) + sd(w.score2.dp$Freq), color="blue", linetype = "longdash") +
  geom_hline(yintercept = mean(w.score2.dp$Freq) - sd(w.score2.dp$Freq), color="blue", linetype = "longdash") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  #ylim(0, 2.3) + 
  ylab("Weighted score") + 
  xlab("Cow ID") + 
  coord_flip()

# Plot w.score 
w.score1.dp <- data.frame(w.score)

w.score1.dp$Var1 <- factor(w.score1.dp$Var1, levels=levels(w.score2.dp$Var1))

w.score1.p <- ggplot(w.score1.dp, aes(x=Var1, y=Freq)) + 
  geom_bar(stat="identity") + 
  theme_bw() + 
  geom_hline(yintercept = mean(w.score1.dp$Freq), color="blue") +
  geom_hline(yintercept = mean(w.score1.dp$Freq) + sd(w.score1.dp$Freq), color="blue", linetype = "longdash") +
  geom_hline(yintercept = mean(w.score1.dp$Freq) - sd(w.score1.dp$Freq), color="blue", linetype = "longdash") + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) + 
  #ylim(0, 2.3) + 
  ylab("Weighted score") + 
  xlab("Cow ID") + 
  coord_flip()

# Combine and save plots
w.score.pcomb <- plot_grid(w.score2.p, w.score1.p)

ggsave(filename = "output/1_shedding_analysis/w_score_comb.pdf", w.score.pcomb, width = 7, height = 3.5)

ggsave(filename = "output/1_shedding_analysis/w_score_comb.jpg", w.score.pcomb, width = 7, height = 3.5, dpi = 400)

ggsave(filename = "output/1_shedding_analysis/w_score_main.pdf", w.score2.p, width = 4, height = 3.5)

ggsave(filename = "output/1_shedding_analysis/w_score_main.jpg", w.score2.p, width = 4, height = 3.5, dpi = 400)
