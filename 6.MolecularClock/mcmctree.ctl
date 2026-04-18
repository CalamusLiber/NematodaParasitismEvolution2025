          seed = -1
       seqfile = seq.phy
      treefile = input.tre
       outfile = out

         ndata = 1
       seqtype = 2  * 0: nucleotides; 1:codons; 2:AAs
       usedata = 3    * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV
         clock =     * 1: global clock; 2: independent rates; 3: correlated rates
       RootAge =   * safe constraint on root age, used if no fossil for root.

         model = 0    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85 (no effect when using the likelihood approximation)
         alpha = 0    * alpha for gamma rates at sites (no effect when using the likelihood approximation)
         ncatG = 4    * No. categories in discrete gamma (no effect when using the likelihood approximation)

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

*** The following have no effect in step 1 when using the likelihood approximation

       BDparas = 1 1 0 M   * birth, death, sampling, conditional/multiplicative, manually change the last number when setting sample fraction as 0.005 or 1.
   kappa_gamma = 6 2      * gamma prior for kappa (no effect when using the likelihood approximation)
   alpha_gamma = 1 1      * gamma prior for alpha (no effect when using the likelihood approximation)

   rgene_gamma = 2 x 1   * gamma prior for overall rates for genes, the x will be replaced by bash script.
  sigma2_gamma = 1 10 1    * gamma prior for sigma^2     (for clock=2 or 3)

      finetune = 1: 0.1  0.1  0.1  0.1  0.1  0.1  * auto (0 or 1) : times, musigma2, rates, mixing, paras, FossilErr (set as auto = 1)

         print = 1
        burnin = 200000
      sampfreq = 15
       nsample = 200000
	
*** Note: Make your window wider (100 columns) before running the program.
