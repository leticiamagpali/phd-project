#!/bin/bash 
#$ -S /bin/bash 
. /etc/profile 
#$ -o /home/leticiamagpali/phd/evol_models/codeml/small_dataset/runs_models_test 
#$ -M leticiamagpali@dal.ca
#$ -m be 
cd /home/leticiamagpali/phd/scripts 
source activate paml
python RunCodeml.py
