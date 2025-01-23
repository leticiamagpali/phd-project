####### Label phylohgenetic trees #########
###########################################

# Author: Letícia Magpali
''' This script labels tree nodes and branches.
To run it, type python label-trees.py <foreground branches> <label> <path to trees folder>
'''

import os
import sys
from re import sub
from ete3 import EvolTree


class EvolTree2(EvolTree):

    def write(self, features=None, outfile=None, format=13):
        if format == 13:
            nwk = super().write(format=11)
            nwk = sub(r'\[&&NHX:mark=([ #\$0-9.]*)\]', r'\1', nwk)
        else:
            nwk = super().write(features=features, format=format)

        if outfile is not None:
            with open(outfile, "w", encoding="UTF-8") as output_file:
                output_file.write(nwk)
            return nwk
        else:
            return nwk

    def get_descendants_by_name(self, name: str):
        for node in self.traverse():
            if node.name == name:
                return node
        return None


# path to file with foreground clades Asin,Psin:Pcat,Dleu:Ddel
FOREGROUND = sys.argv[1]
LABEL = sys.argv[2]  # label
HYPOTHESIS = sys.argv[3]
TREES_FOLDER = sys.argv[4]  # insert here path to trees folder
OUTPUT_FOLDER = sys.argv[5]

# Takes user input, a string of with the names of foreground branches
# separated by spaces, and turns it into a list
foreground_nodes = FOREGROUND.split(",")

# This loop goes through each tree in your trees folder,
# parses the tree and checks only for the leafs
# if a particular leaf is in the foreground list,
# it will add the user specified label
for file in os.listdir(TREES_FOLDER):
    if file.endswith(".tre"):
        gene_name = file.split("_tree_M0.tre")[0]
        tree_path = os.path.join(TREES_FOLDER, file)
        tree = EvolTree2(tree_path)

        for element in foreground_nodes:
            if ":" not in element:
                # this is a leaf
                leaf = tree.get_descendants_by_name(element)
                if leaf is not None:
                    tree.mark_tree([leaf.node_id], marks=[f"#{LABEL}"])
            else:
                # this is a clade
                clade_limits = element.split(":")
                lim_left = tree.get_descendants_by_name(clade_limits[0])
                lim_right = tree.get_descendants_by_name(clade_limits[1])
                common = lim_left.get_common_ancestor(lim_right)
                tree.mark_tree([common.node_id], marks=[f"${LABEL}"])

        tree.write(outfile=f"{OUTPUT_FOLDER}/{gene_name}_tree_{HYPOTHESIS}.tre", format=13)
