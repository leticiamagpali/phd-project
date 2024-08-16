#!/bin/bash 
#$ -S /bin/bash 
. /etc/profile 
#$ -o /home/leticiamagpali/evol_models/codeml/branch_models/M0/M0_v3
#$ -M leticiamagpali@dal.ca
#$ -pe threaded 10 
cd /home/leticiamagpali/evol_models/codeml/branch_models/M0/M0_v3
source activate paml
codeml M0-v3.ctl
