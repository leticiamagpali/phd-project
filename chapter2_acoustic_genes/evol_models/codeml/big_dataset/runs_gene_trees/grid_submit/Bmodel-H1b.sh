#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -o /home/leticiamagpali/phd/evol_models/codeml/big_dataset/runs_gene_trees/grid_out
#$ -e /home/leticiamagpali/phd/evol_models/codeml/big_dataset/runs_gene_trees/grid_out
#$ -M leticiamagpali@dal.ca
#$ -m be
cd /home/leticiamagpali/phd/evol_models/codeml/big_dataset
source activate paml
python run-codeml.py /home/leticiamagpali/phd/evol_models/codeml/big_dataset/runs_gene_trees/Bmodel-H1b
