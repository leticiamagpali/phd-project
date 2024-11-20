#!/bin/bash 
#$ -S /bin/bash 
. /etc/profile 
#$ -o /home/leticiamagpali/phd/evol_models/codeml/small_dataset/runs_models
#$ -e /home/leticiamagpali/phd/evol_models/codeml/small_dataset/runs_models
#$ -M leticiamagpali@dal.ca 
#$ -m be 
cd /home/leticiamagpali/phd/evol_models/codeml/small_dataset/runs_models/Amodel-H3
source activate paml
python RunCodeml_Amodel-H3.py
