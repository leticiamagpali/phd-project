      seqfile = ../../../cml_align/codon-align-h1-cml.phy                                  * Path to the alignment file
     treefile = ../../../cml_trees/speciestrees-h1b.phy                                        * Path to the tree file
      outfile = out_Anull_h1b.txt                                                                * Path to the output file
   
        noisy = 3                     * How much rubbish on the screen
      verbose = 1                     * More or less detailed report

      seqtype = 1                      * Data type
        ndata = 9 separate_trees       * Number of data sets or loci
        icode = 0                      * Genetic code 
    cleandata = 0                      * Remove sites with ambiguity data?
		
        model = 2                     * Models for ω varying across lineages
	  NSsites = 2                       * Models for ω varying across sites
    CodonFreq = 7                     * Codon frequencies
	  estFreq = 0                       * Use observed freqs or estimate freqs by ML
        clock = 0                     * Clock model
    fix_omega = 1                     * Estimate or fix omega
        omega = 1                     * Initial or fixed omega