#!bin/bash

# The following script will run all anlysis. 
# This will take a long time. 

yourdir=$(pwd)

source activate $HOME/local_conda_env/r_clone/

Rscript $yourdir/scr/1_shedding_analysis_conv.R 

Rscript 

conda deactivate 