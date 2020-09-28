#!bin/bash 

source activate /home/umane001/local_conda_env/r_3.6.1/

Rscript /home/umane001/WBVR_MAP_Microbiota/scr/rfPremute_run.R

conda deactivate 