#This sript will create all nessesary folder for a tiypical R project. 
#Also script will create a file for github to ignor all files but scr folder 
#Prior runing the script set working directory (setwd()) to your R project folder. 
#The script will add only not extisting directories 
####################################################

#1. Add directories 
l.dir <- c("scr", "data", "output", "resourses", 
           "scr/functions", 
           "data/Dada2_input", "data/metadata", 
           "output/plots", "output/objects")

for (d in l.dir) {
    dir.create(d, showWarnings = FALSE)
}

#2. Add .gitignore file with. 
#   All folders but scr are ignored by defult 
fileConn <- file("gitignore.txt")
writeLines(c("data", "output", "resourses", 
           "data/Dada2_input", "data/metadata", 
           "output/plots", "output/objects", 
            ".ipynb_checkpoints"), fileConn)
close(fileConn)

file.rename("gitignore.txt", ".gitignore")