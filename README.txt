# The following analysis pipeline is supplementary material for the publication: “Faecal microbiota composition based random forest model predicts Mycobacterium Avium subsp. Paratuberculosis (MAP) shedding severity in cattle”
 
 The goal of this pipeline is to provide a clear insight into statistical analysis and
 promote reproducibility of results.
 
 Things to keep in mind:
 1. The pipeline is not intended as a guideline for analysis of microbiota related data and
 contains a large number of functions applicable only in the context of the data used for publications.

 2. The pipeline is not optimized for performance on weak machines and in particular for RAM use.Certain scripts will require around considerable (32  Gb or more) of RAM to run and could crash if a system has less. 

 3. To follow the pipeline step by step a Linux operating system (Debain) is required and if another system is used it should be addressed accordingly.
 

 # Analysis reproduction #
 ####################

# 1. Clone the gitHub repository (please consult a tutorial regarding installation and use of gitHub tools).

git clone https://github.com/AlexanderUm/WBVR_MAP_Microbiota.git


# 2. Navigate to the cloned folder and run script to make directories

cd WBVR_MAP_Microbiota/

Rscript make_project_structure.R


# 3. Copy fastq.gz files into data/Dada2_input
Data is available from NCBI SRA archive: 


# 4. Copy sedding data into data/
Shedding data is available as supplementary material. 
File with shedding data should be named “shedding_data.csv”


# 5. Copy metadata into data/metadata
metadata data is available as supplementary material. 
File with metadata should be named “samples_metadata_f.csv”


# 6. Copy reference database (see DADA2 tutorial) for taxonomy assignment to resources/ folder
We used silva_nr_v138_train_set.fa.gz

# 7. Scripts are ready to run. 
You can run them from jupyter notebooks (files with the extension .ipynb) or using R scripts (files with the extension .R). 
Files should be run consequently 1_xxx, 2_xxx, 3_xxx, ….

