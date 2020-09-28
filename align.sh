#!/bin/bash 
for i in $HOME/WBVR_MAP_Microbiota/output/fasta_for_align/*
do
muscle -in $i -out $i
done
 