      seqfile = ../../../../cml_align/codon-align-v4-cml.phy                                    * Path to the alignment file
     treefile = ../../../../cml_trees/speciestree-v4.phy                                        * Path to the tree file
      outfile = out_M0_v4.txt                                                                * Path to the output file
   
        noisy = 3                     * How much rubbish on the screen
      verbose = 1                     * More or less detailed report

      seqtype = 1                      * Data type
        ndata = 10 maintree            * Number of data sets or loci
        icode = 0                      * Genetic code 
    cleandata = 0                      * Remove sites with ambiguity data?
		
        model = 0                     * Models for ω varying across lineages
	  NSsites = 0                       * Models for ω varying across sites
    CodonFreq = 3                     * Codon frequencies
	  estFreq = 0                       * Use observed freqs or estimate freqs by ML
        clock = 0                     * Clock model
    fix_omega = 0                     * Estimate or fix omega
        omega = 0.5                   * Initial or fixed omega