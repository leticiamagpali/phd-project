###################################
# Generates control file for codeml
###################################

# Author: Letícia Magpali

'''This script generates a control file for codeml by editing a 
file template and making multiple copies of it for each gene.
Each control file for each gene/analysis is stored in its own folder.
'''

import os
import sys
from pathlib import Path

# Text file with analysis name and parameters -> taken as input
# Control file template taken as a input
# Hypothesis taken as a input

# Insert here your file locations for codeml

RUN_FOLDER = sys.argv[1] #path to runfolder
PATHS_FILE = sys.argv[2] # file with paths to alignments folder and tree folder one per line
PARAMETERS_FILE = sys.argv[3] # file with analysis name and a list of codeml parameters, one per line
HYPOTHESIS = sys.argv[4]  # number/code of your hypothesis to be written in folders and outfiles
CONTROL_FILE_TEMPLATE = sys.argv[5] # path or name of control file template

# Reading input files and assigning each line of the file to a path variable
with open(PATHS_FILE, "r") as paths_file:
    paths = [line.strip() for line in paths_file]

ALIGNMENTS_FOLDER, TREES_FOLDER = paths[:2]

# Converting to path
RUN_FOLDER = Path(RUN_FOLDER)
ALIGNMENTS_FOLDER = Path(ALIGNMENTS_FOLDER)
TREES_FOLDER = Path(TREES_FOLDER)

# Print to verify
print(f"run folder: {RUN_FOLDER}")
print(f"alignments folder: {ALIGNMENTS_FOLDER}")
print(f"treees folder: {TREES_FOLDER}")

with open(PARAMETERS_FILE, "r") as parameters_file:
    parameters = [line.strip() for line in parameters_file]

ANALYSIS, BRANCH_MODEL, SITE_MODEL, CODON_FREQ_MODEL, OMEGA, FIX_OMEGA = parameters[:6]

NCATG = None

if len(parameters) > 6: 
    NCATG = parameters[6]
else:
    None 

print(f"model = {BRANCH_MODEL}")
print(f"NSsites = {SITE_MODEL}")
print(f"CodonFreq = {CODON_FREQ_MODEL}")
print(f"omega = {OMEGA}")
print(f"fix_omega = {FIX_OMEGA}")
print(f"ncatG = {NCATG}")

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
        
        # Add `NCAT` if it’s not None
        if NCATG is not None:
            codeml_dictionary["ncatG"] = NCATG

        # Creates a folder to store the ctl file
        run_subfolder = RUN_FOLDER/f"{ANALYSIS}-{HYPOTHESIS}_{gene_name}"

        run_subfolder.mkdir(parents=True, exist_ok=True)

        # and writes the dictionary to a text file to create a control file
        with open(run_subfolder/f"{ANALYSIS}-{HYPOTHESIS}_{gene_name}.ctl", "w", encoding="utf-8") as control_file:
            for key, val in codeml_dictionary.items():
                control_file.write(f"{key} = {val}\n")
