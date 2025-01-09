#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -o /home/leticiamagpali/evol_models/codeml/branch_models/h1_whales_dolphins/2M_h1a
#$ -e /home/leticiamagpali/evol_models/codeml/branch_models/h1_whales_dolphins/2M_h1a
#$ -M leticiamagpali@dal.ca
#$ -m be
#$ -pe threaded 4
cd /home/leticiamagpali/evol_models/codeml/branch_models/h1_whales_dolphins/2M_h1a
source activate paml
codeml h1a.ctl
