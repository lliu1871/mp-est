# MP-EST (v2.1): Maximum Pseudo-likelihood Estimation of Species Trees
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

MP-EST estimates species trees (topology and branch lengths in coalescent units) from a set of gene trees by maximizing a pseudo-likelihood function. The input data of MP-EST are rooted gene trees. In addition to the gene tree file, a control file must be generated for running MP-EST. The control file contains necessary parameters for running MP-EST.

## New features in version 2.1

1. The program can run on millions of rooted gene trees (very fast). The computational time does not rely on the number of gene trees.

2. The program can take polytomy gene trees as input to estimate species trees.

3. The program can convert short branches (<1e-06) to polytomies in the gene trees.

4. The program can calculate triplet distances among gene trees

5. The program can calculate quartet distances among gene trees

6. The program outputs the species tree gene tree concordance scores for each internal node of the MP-EST tree.


## Compile from source code
To compile the program from source code, type make and hit return under the directory src.

## Run the program
./mpest controlfile

## Example control files
Example control files are included in the example folder. 

A control file for building the MPEST tree (control)

testgenetree #gene tree file

0   #0:running mp-est, 1: calculating triplet distances, 2: converting short branches to polytomies, 3: calculating quartet distances, 4: calculating gene tree partitions

-1  # seed, -1: a random seed

2   # the number of runs

100 5 #the number of genes and the number of species

A 1 A #species-allele table

B 1 B

C 1 C

D 1 D

E 1 E

0 #0:no user tree, 1: the user tree is used as the start tree in the algorithm

(E:0.002217212018,((C:0.001601580862,(A:0.0007794028167,B:0.001370162343):0.0004631267076):0.0003658417881,D:0.001319072462):0.0003568793102); #user tree


## Ouput files
There are two output files; testgenetree_besttree.tre and testgenetree_output.tre. The trees updated in the algorithm are saved in testgenetree_output.tre, while the mpest tree is saved in testgenetree_besttree.tre. If multiple runs are specified in the control file, testgenetree_besttree.tre contains multiple mpest trees.


## Citation
Liu, L., L. Yu, S.V. Edwards. A maximum pseudo-likelihood approach for estimating species trees under the coalescent model. BMC Evol. Biol. 2010, 10:302.


## Old versions
Old versions MP-EST are available at https://github.com/lliu1871/oldversion
