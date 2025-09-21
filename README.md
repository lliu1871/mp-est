# MP-EST (v3.0): Maximum Pseudo-likelihood Estimation of Species Trees

[![Install](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://bioconda.github.io/recipes/mp-est/README.html)
[![OS](https://anaconda.org/bioconda/mp-est/badges/platforms.svg)](#mp-est)
[![Version](https://img.shields.io/conda/vn/bioconda/mp-est?label=version)](https://bioconda.github.io/recipes/mp-est/README.html)
[![Release Date](https://anaconda.org/bioconda/mp-est/badges/latest_release_date.svg)](#mp-est)
[![Downloads](https://img.shields.io/conda/dn/bioconda/mp-est.svg?style=flat)](https://bioconda.github.io/recipes/mp-est/README.html)
[![License](https://anaconda.org/bioconda/mp-est/badges/license.svg)](https://github.com/lliu1871/mp-est/blob/master/LICENSE)

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

MP-EST estimates species trees (topology and branch lengths in coalescent units) from a set of gene trees by maximizing a pseudo-likelihood function. The input data of MP-EST are rooted gene trees. Unlike previous versions, MP-EST3.0 uses the command line interface, i.e., no control file is needed for running MP-EST3.0. 

## New features in version 3.0

1. MP-EST3.0 can build an NJst tree and use it as the initial tree to find the MP-EST tree.

2. MP-EST3.0 can convert short branches (<1e-06) to polytomies in the gene trees.

3. MP-EST3.0 can take polytomy gene trees as input to estimate species trees.

4. MP-EST3.0 can calculate triplet/quartet distances among gene trees

5. MP-EST3.0 outputs species-tree-gene-tree concordance scores for each internal node of the MP-EST tree.

# Installation

There are two ways to install MP-EST:

## 1. Install from [bioconda](https://bioconda.github.io/recipes/mp-est/README.html)

MP-EST is available as a package on bioconda. With [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) or [mamba](https://mamba.readthedocs.io/en/latest/installation.html) already installed and your [channels set up for bioconda](https://bioconda.github.io/#usage), simply type the following to install:

```bash
conda create -n mp-est-env
conda activate mp-est-env
conda install mp-est
```

This will create a new environment and activate it and then install mp-est from bioconda. Alternateively, you could install in an environment you've already created.

If you use mamba instead of conda, just replace `conda` with `mamba` above.

## 2. Compile from source code

To compile the program from source code, type make and hit return under the directory src.

# Usage

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
