"""
This script removes short isomorphs through the keep_longest_isoform function which
receives a filepath to some FASTA file containing a collection of sequences from some
species. These sequences are read into a list which is sorted by descending length.
The headers of the sequence are then checked for a distinguishing attribute value, 
which identifies to which gene/protein/locus_tag that sequence belongs.
Novel values are added to a running list, and the associated sequences are appended
to make the new dataset which is returned at the end. This assures that:
1) Each copy of a gene (represented by the attribute value) only appears once in the new dataset
2) The sequence is the longest representative for that attribute

To add compatibility for more species, expand the dictionary in the 
get_distinguishing_attribute function

To utilize this tool, access through run.py

The distinguishing attribute is used to identify each gene/protein. 
If that gene/protein has already been included in the output, that means
that the longest isoform has aleady been extracted and other copies of the
gene will not be included. 
"""

from Bio import SeqIO
import os.path

# Return the list of sequence records in descending order


def sortSeq_descending_size(sequences):
    return sorted(sequences, key=lambda r: -len(r))


######## Identifying gene tags ########

'''This function uses a dictionary that links the name of the folder containing each genome 
to the gene tags used in that genome's annotation. 
Modify the dictionary so that: 
    - keys = your folder names, 
    - values = list of all gene tags (this may vary depending on how each genome was annotated)
'''


def get_distinguishing_attribute(token):
    if os.path.exists(str(token)):
        filepath = token
        # To add compatibility for more species or gene tags, expand this dictionary
        # "gene=", "locus_tag="
        label = {}  # how to make dictionary from a list
        return str(label[filepath.split('/')[-3]])
    else:
        seq_record = token
        if "isoform".casefold() in seq_record.description.split(" ")[1:-2]:
            label = " ".join(seq_record.description.split(" ")[1:-3])
        else:
            label = " ".join(seq_record.description.split(" ")[1:-2])
        return str(label)

# Collect all sequences annotated to the same gene or protein


def keep_longest_isoform(filepath):
    attributes_added = []
    new_dataset = []
    # Create a list of sequences sorted by descending length
    sequences = sortSeq_descending_size(list(SeqIO.parse(filepath, "fasta")))
    # Determine the distnguishing attribute name for nucleotide files
    if ".fna" in filepath:
        distinguishing_attribute = get_distinguishing_attribute(filepath)
        # Iterate through records adding the associated record of each novel enounter of an attribute to the new dataset
        for seq_record in sequences:
            attribute_value = [s for s in seq_record.description.split(
                " ") if distinguishing_attribute in s]
            if not (attribute_value in attributes_added):
                attributes_added.append(attribute_value)
                new_dataset.append(seq_record)
        return new_dataset
    # Determine the distnguishing attribute name for protein files
    else:
        # Iterate through records adding the associated record of each novel enounter of an attribute to the new dataset
        for seq_record in sequences:
            attribute_value = get_distinguishing_attribute(seq_record)
            if not (attribute_value in attributes_added):
                attributes_added.append(attribute_value)
                new_dataset.append(seq_record)
        return new_dataset


# makes sure functions will work when called from this file AND also when called from run.py
if __name__ == '__main__':
    keep_longest_isoform()


# This code was adapted from code originally written by Yuri on [date]
# written collaborativelly for [project name]
