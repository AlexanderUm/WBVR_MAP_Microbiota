#!bin/bash

# The following script will run all anlysis. 
# This will take a long time. 

yourdir=$(pwd)

source activate $HOME/local_conda_env/r_clone/

# Script 1 
Rscript $yourdir/scr/1_shedding_analysis_conv.R 

echo script 1 complitted > output/log.txt

# Script 2 
Rscript $yourdir/scr/2_dada2_processing.R

echo script 2 complitted >> output/log.txt

# Script 3
Rscript $yourdir/scr/3_filtering_normalization.R 

echo script 3 complitted >> output/log.txt

# Script 4 
Rscript $yourdir/scr/4_betta_diversity.R 

echo script 4 complitted >> output/log.txt

# Script 5
Rscript $yourdir/scr/5_RF_regression_model.R 

echo script 5 complitted >> output/log.txt

# Script 6
Rscript $yourdir/scr/6_RF_model_general_best.R 

echo script 6 complitted >> output/log.txt

# Script 7 
Rscript $yourdir/scr/7_discriminatory_tax.R

echo script 7 complitted >> output/log.txt

# Script 8 
Rscript $yourdir/scr/8_lean_RF_model.R 

echo script 8 complitted >> output/log.txt

Rscript $yourdir/scr/9_time_influence/9_time_influence_1.R

Rscript $yourdir/scr/9_time_influence/9_time_influence_2.R

Rscript $yourdir/scr/9_time_influence/9_time_influence_3.R

Rscript $yourdir/scr/9_time_influence/9_time_influence_4.R

Rscript $yourdir/scr/9_time_influence/9_time_influence_5.R

Rscript $yourdir/scr/9_time_influence/9_time_influence_6.R

Rscript $yourdir/scr/9_time_influence/9_time_influence_7.R

Rscript $yourdir/scr/9_time_influence/9_time_influence_8.R

Rscript $yourdir/scr/9_time_influence/9_time_influence_9.R

Rscript $yourdir/scr/9_time_influence/9_time_influence_10.R

Rscript $yourdir/scr/9_time_influence/9_time_influence_11.R

Rscript $yourdir/scr/9_time_influence/9_time_influence_12.R


conda deactivate 