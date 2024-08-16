###################################
# Generates control file for codeml
###################################

# Author: Letícia Magpali

'''This script generates a control file for codeml by editing a 
file template and making multiple copies of it for each gene.
Each control file for each gene/analysis is stored in its own folder.
'''

import os
from pathlib import Path

# Insert here your file locations for codeml

# Edit this one time, and it will stay the same as long as you are working with the same dataset
ALIGNMENTS_FOLDER = Path(
    "/home/leticiamagpali/phd/evol_models/codeml/small_dataset/cml_align")

# Change this every time you are running a different hypothesis (every hypothesis has one different LABELLED tree)
TREES_FOLDER = Path(
    "/home/leticiamagpali/phd/evol_models/codeml/small_dataset/cml_trees/H2_trees")

# Edit this one time only, to tell the script where your control file template is
CONTROL_FILE_TEMPLATE = Path(
    "/home/leticiamagpali/phd/evol_models/codeml/small_dataset/control_file_cml.ctl")

# Change this for every different run
# this is the folder where all your runs for different genes are, for that particular model
RUN_FOLDER = Path(
    "/home/leticiamagpali/phd/evol_models/codeml/small_dataset/runs_models/Amodel-null-H2")

# Change this with the name of the model you are running
ANALYSIS = "2model"

# Change this with the number/code of your hypothesis
HYPOTHESIS = "H2"


# Insert here your parameters for codeml, according to your model
BRANCH_MODEL = "2"
SITE_MODEL = "0"
CODON_FREQ_MODEL = "7"
OMEGA = "1"
FIX_OMEGA = "0"

######## RUNNING CODE - NOT FOR EDITING ##############
######################################################

# Creates a dictionary of the codeml template file
codeml_dictionary = {}
with open(CONTROL_FILE_TEMPLATE) as control_file_template:
    for line in control_file_template:
        parametro, comentario = line.split("*")
        (key, val) = parametro.split("=")
        key = key.strip()
        val = val.strip()
        codeml_dictionary[key] = val

# This loop creates a control file from the dictionary
# and stores it on a run folder
for alignment in os.listdir(ALIGNMENTS_FOLDER):
    if alignment.endswith(".phy"):
        alignment_path = ALIGNMENTS_FOLDER/alignment
        # Extracts gene name from alignment file name
        gene_name = alignment.split("_codon_aligned.phy")[0]
        # Take tree path from trees folder
        tree_path = TREES_FOLDER/f"{gene_name}_tree_{HYPOTHESIS}.tre"

        # Modifies dictionary adding the values you specified above
        codeml_dictionary["seqfile"] = alignment_path
        codeml_dictionary["treefile"] = tree_path
        codeml_dictionary["outfile"] = f"out_{gene_name}_{ANALYSIS}-{HYPOTHESIS}.txt"
        codeml_dictionary["omega"] = OMEGA
        codeml_dictionary["fix_omega"] = FIX_OMEGA
        codeml_dictionary["model"] = BRANCH_MODEL
        codeml_dictionary["NSsites"] = SITE_MODEL
        codeml_dictionary["CodonFreq"] = CODON_FREQ_MODEL

        # Creates a folder to store the ctl file
        run_subfolder = RUN_FOLDER/f"{ANALYSIS}-{HYPOTHESIS}_{gene_name}"

        run_subfolder.mkdir(parents=True, exist_ok=True)

        # and writes the dictionary to a text file to create a control file
        with open(run_subfolder/f"{ANALYSIS}-{HYPOTHESIS}_{gene_name}.ctl", "w", encoding="utf-8") as control_file:
            for key, val in codeml_dictionary.items():
                control_file.write(f"{key} = {val}\n")
