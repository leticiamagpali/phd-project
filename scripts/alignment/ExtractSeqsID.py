#### Editing alignments ####
############################

'''This module has functions to edit multispecies alignments'''

# Importing necessary packages
from Bio import SeqIO
from ..file_editing import TextToList


def extract_seqs_id(alignment_path, list_path):
    '''This function extracts sequences from a phylip file
    based on a list of the sequence IDs'''

    # Creates a list based on the sequence
    seqs_list = TextToList.create_list(list_path)

    # Transforms seqs_list into a set
    seqs_set = set(seqs_list)

    # creates empty list to store the sequence records of the extracted sequences
    # extracted_seqs = MultipleSeqAlignment([])
    extracted_seqs = []

    # reads sequence file
    for record in SeqIO.parse(open(alignment_path), "fasta"):

        # if sequence ID matches the name on seqs_list
        # appends sequence record to the extracted_seqs list
        if record.id in seqs_set:
            extracted_seqs.append(record)
    return (extracted_seqs)


# Make changes to make it work on both alignments and sequences
