for i in $(find codon_alignments_edited/ -name "*_codon_aligned.phy"); do iqtree -s ${i} -st CODON -m TEST -o Mmus_${i:24:(${#i}-24-18)} -pre ${i:24:(${#i}-24-18)} -nt 4; done
