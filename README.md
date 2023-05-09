# MP-EST (v3.0): Maximum Pseudo-likelihood Estimation of Species Trees
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

MP-EST estimates species trees (topology and branch lengths in coalescent units) from a set of gene trees by maximizing a pseudo-likelihood function. The input data of MP-EST are rooted gene trees. In addition to the gene tree file, a control file must be generated for running MP-EST. The control file contains necessary parameters for running MP-EST.

## New features in version 3.0

1. The version 3.0 uses the command line interface. No control file is needed.

2. The version 3.0 can build an NJst tree and use it as the initial tree to find the MP-EST tree.

3. The program can convert short branches (<1e-06) to polytomies in the gene trees.

4. The program can take polytomy gene trees as input to estimate species trees.

5. The program can calculate triplet/quartet distances among gene trees

6. The program outputs the species tree gene tree concordance scores for each internal node of the MP-EST tree.


## Compile from source code
To compile the program from source code, type make and hit return under the directory src.

## Run the program
./mpest -i testgenetree 

## Help
./mpest -h

Usage: mpest [-i inputfile] [-n #] [-s #] [-u NAME] [-h] [-B|-C|-L|-N|-P|-Q|-T|] \
  -i: NAME = name of input file in nexus format. Input gene trees must be rooted; polytomy trees are allowed.\
  -n: # = number of runs [default = 1].\
  -s: # = seed for random number generator [default = system generated].\
  -u: NAME = name of user tree file.\
  -c: # = convert short branches (<#) to polytomies in gene trees.\
  -h: help message [default = random].\
  -B: optimize branch lengths of a fixed species tree provided through usertree option -u.\
  -L: calculate loglikelihood of a species tree provided through usertree option -u.\
  -N: build NJst tree.\
  -P: calculate partitions for gene trees.\
  -Q: calculate pairwise quartet distances among gene trees. Polytomy (unrooted or rooted) gene trees are allowed.\
  -T: calculate pairwise triple distances of gene trees. Polytomy rooted gene trees are allowed.\

## Ouput files
There are two output files; testgenetree_besttree.tre and testgenetree_output.tre. The trees updated in the algorithm are saved in testgenetree_output.tre, while the mpest tree is saved in testgenetree_besttree.tre. If multiple runs are specified in the control file, testgenetree_besttree.tre contains multiple mpest trees.


## Citation
Liu, L., L. Yu, S.V. Edwards. A maximum pseudo-likelihood approach for estimating species trees under the coalescent model. BMC Evol. Biol. 2010, 10:302.


## Old versions
Old versions MP-EST are available at https://github.com/lliu1871/oldversion
