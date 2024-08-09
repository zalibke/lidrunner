#!/bin/bash
#PBS -l select=1:ncpus=64:mem=100gb
#PBS -l walltime=20:00:00
#Zane Libke
#10-July-2024
#lidrunner v1.4
#this version will keep searching for windows which fit the 
#sites_per_window threshold until LDhat has ran the number specified by windows_per_chunk

#source $HOME/anaconda3/bin/activate

#cd $HOME/lidrunner_coluzzii_0710/code

#Rscript 0-install_rpackages.r

python3 1-runall.py 100 20 500 10 50 TRUE 50 "(taxon == 'coluzzii') and (sex_call == 'F')"





