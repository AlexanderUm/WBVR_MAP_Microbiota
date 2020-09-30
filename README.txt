 # The following analysis pipeline is a suplementary material for the publication. 
 
 The goal of this pipeline is to provide a clear inside into statistical analysis and 
 promote reproducibility of results. 
 
 Disclaimer: 
 1. The pipeline is not intended as a guidline for analysis of microbiota data and 
     contains a large number of fuctions apliciable only in the context in of the used data. 
 2. The pipeline is not optimized for perforamce and in particullary for RAM use. 
     Certain scripts will require around 126 Gb of RAM to run and will crash if a system has less. 
 3. To follow the pipeline step by step a Linux operating system (Debain) required. 
 
 ########################
 # Analysis reproduciton# 
 ########################

# 1. Clone the gitHub repsitory (please consult a tutorial reguarding installation and use of gitHub tools). 

git clone https://github.com/AlexanderUm/WBVR_MAP_Microbiota.git

# 2. Run script to make directories 

cd WBVR_MAP_Microbiota/ 

Rscript make_project_structure.R 

# 3. Copy fastq.gz files into data/Dada2_input

!!!! Destcribe how to pull it from the database !!!

# 4. Copy sedding data into data/

# 5. Copy metadata into data/metadata

# 6. Copy refference database for taxonamy assignment to resourses/ folder 
#    We are using silva_nr_v138_train_set.fa.gz 

# 6. Make run_all.sh executable 

chmod +x run_all.sh 

# 7. Run the script 

bash run_all.sh