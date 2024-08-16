#!/bin/bash 
#$ -S /bin/bash 
. /etc/profile 
#$ -o /home/leticiamagpali/evol_models/codeml/branch_models/M0/M0_v4
#$ -M leticiamagpali@dal.ca
#$ -pe threaded 10
cd /home/leticiamagpali/evol_models/codeml/branch_models/M0/M0_v4
source activate paml
codeml M0-v4.ctl
