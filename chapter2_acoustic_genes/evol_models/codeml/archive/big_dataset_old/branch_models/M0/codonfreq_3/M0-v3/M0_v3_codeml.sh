#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -o /home/leticiamagpali/evol_models/codeml/branch_models/M0/codonfreq_3/M0-v3
#$ -e /home/leticiamagpali/evol_models/codeml/branch_models/M0/codonfreq_3/M0-v3
#$ -M leticiamagpali@dal.ca
#$ -m be
#$ -pe threaded 4
cd /home/leticiamagpali/evol_models/codeml/branch_models/M0/codonfreq_3/M0-v3
source activate paml
codeml M0-v3.ctl
