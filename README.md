# MP-EST (v2.1): Maximum Pseudo-likelihood Estimation of Species Trees
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

MP-EST estimates species trees from a set of gene trees by maximizing a pseudo-likelihood function. The input data of MP-EST are rooted binary gene trees produced by the maximum likelihood phylogenetic programs RAxML, PHYML, PHYLIP, and PAUP etc. In addition to the gene tree file, a control file must be generated for running MP-EST. The control file contains necessary parameters for running MP-EST.

##New features in version 2.1
Parallel computing using PTHREAD is available in version 2.0. Set PTHREAD ?= yes in Makefile
The constraints on the birth rate parameter lambda in the old versions are removed in version 2.0
The birth-death probability is calculated using two formulas - one for alpha < 0.5 and another for alpha > 0.5 where alpha = lambda*brlens / (1+lambda*brlens)
The program provides a warning message if the species names in the data file do not match the names in the control file
The program provides a warning message if the number of gene families in the control file > the number of gene families in the data file
The second column "NUMBER" is removed in the data file
Compile from source code
To compile the program from source code, type make and hit return under the directory src.

Run the program
./begfe controlfile

Example control files
Three control files are included in the package. The control file controlsim is used to simulate gene family data.

A control file for simulating gene family data (controlsim)
1 #0:analysis, 1:simulation

sim1 #output file

-1 #random seed

100 #number of gene families to simulate

5 #number of species

(((chimp:6#0.005,human:6#0.005):81#0.005,(mouse:17#0.005,rat:17#0.005):70#0.005):6#0.005,dog:93#0.005)#0.005; #species tree, the numbers after '#' are birth/death rates

5 15 #the number of gene copies at the tree root is simulated from an uniform distribution (5, 15)

To run the simulation, type ./begfe controlsim

This will produce two files; sim1 and sim1.true. The simulated gene family dataset is saved in the file sim1, while sim1.true contains the true values of the parameters in the Bayesian model.

A control file for analyzing gene family data (control)
The other control file control is for carrying out the Bayesian analysis of the gene family data. Type ./begfe control and hit return.

0 #0:analysis, 1:simulation

sim1 #input file

-1 #random seed

100 #number of gene families

5 #number of species

(((chimp:6,human:6):81,(mouse:17,rat:17):70):6,dog:93); #species tree

10000 100 0 #number of MCMC generations, save every 100 samples, 0:unlinked (variable) lambdas and 1:linked (single) lambda

A second example of analyzing gene family data (AnolisMHC_control.txt)
0 #0:analysis, 1:simulation

AnolisMHCtab.txt #input file

-1 #random seed

4 #number of gene families

13 #number of species

(((((((Laticauda_laticaudata:32.71005600,(Pseudonaja_textilis:27.60000000,Notechis_scutatus:27.60000000):5.11005600):2.28994400,Naja_naja:35.00000000):132.12419034,(Anolis_carolinensis:66.08013667,Anolis_sagrei:66.08013667):101.04405368):10.87262926,(Podarcis_muralis:148.65977538,Salvator_merianae:148.65977538):29.33704422):73.83363773,Sphenodon_punctatus:251.83045733):27.82651933,(Gallus_gallus:98.04286929,Taeniopygia_guttata:98.04286929):181.61410737):32.24694470,(Mus_musculus:89.82318742,Homo_sapiens:89.82318742):222.08073394);

1000000 1000 0 #number of MCMC generations, save every 100 samples, 0:unlinked (variable) lambdas and 1:linked (single) lambda

Ouput files
There are two output files; sim1.out and sim1.pvalue. The MCMC output for all parameters is saved in sim1.out. Each column in sim.out represents the posterior distribution of a particular parameter in the birth and death model. The parameters such as the birth and death rate are estimated by the Bayesian means, i.e., the averages of the columns after discarding the burn-in period.

The Bayesian p-values (PPP) for all gene families in the dataset are saved in sim1.pvalue. The average of a column in sim1.pvalue is the Bayesian p-values for a particular gene family. Of course, the burn-in period must be discarded.

Citation
Liu, L., L. Yu, S.V. Edwards. A maximum pseudo-likelihood approach for estimating species trees under the coalescent model. BMC Evol. Biol. 2010, 10:302.


Old versions
Old versions BEGFE are available at https://github.com/lliu1871/oldversion
