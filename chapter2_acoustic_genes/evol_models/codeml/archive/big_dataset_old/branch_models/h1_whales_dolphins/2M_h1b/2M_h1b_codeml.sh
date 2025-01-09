#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -o /home/leticiamagpali/evol_models/codeml/branch_models/h1_whales_dolphins/2M_h1b
#$ -e /home/leticiamagpali/evol_models/codeml/branch_models/h1_whales_dolphins/2M_h1b
#$ -M leticiamagpali@dal.ca
#$ -m be
#$ -pe threaded 4
cd /home/leticiamagpali/evol_models/codeml/branch_models/h1_whales_dolphins/2M_h1b
source activate paml
codeml h1b.ctl
