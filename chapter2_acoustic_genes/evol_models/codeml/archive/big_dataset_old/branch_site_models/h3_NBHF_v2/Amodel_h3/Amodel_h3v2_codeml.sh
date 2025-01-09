#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -o /home/leticiamagpali/evol_models/codeml/branch_site_models/h3_NBHF_v2/Amodel_h3
#$ -e /home/leticiamagpali/evol_models/codeml/branch_site_models/h3_NBHF_v2/Amodel_h3
#$ -M leticiamagpali@dal.ca
#$ -m be
#$ -pe threaded 4
cd /home/leticiamagpali/evol_models/codeml/branch_site_models/h3_NBHF_v2/Amodel_h3
source activate paml
codeml A_model_h3.ctl
