source activate iqtree

#ModelFinder
for i in $(find codon_alignments_v*_Jan24/ -name "*_codon_aligned.phy"); do iqtree -s ${i} -st CODON -m TEST -pre ${i:26:(${#i}-26-18)}_MF -nt 4; done

#GTR + G4
for i in $(find codon_alignments_v*_Jan24/ -name "*_codon_aligned.phy"); do iqtree -s ${i} -m GTR+G4 -pre ${i:26:(${#i}-26-18)}_GTR -nt 4; done

#WAG + G4
for i in $(find translated/ -name "*_codon_aligned.translated.fasta"); do iqtree -s ${i} -st AA -m WAG+G4 -pre ${i:11:(${#i}-11-31)}_WAG -nt 4; done
