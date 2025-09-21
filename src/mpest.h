#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <time.h>
#include <stdarg.h>
#include <sys/time.h>
#include <unistd.h>

#define MAX_STRING_LENGTH 10000 /*max length of a string*/
#define NTAXA         200      /* max # of species */
#define MAXROUND	10000000		/* MAX # OF ROUNDS*/
#define NUM_NOCHANGE	20000		/* # OF ROUNDS THAT NO BIGGER LIKELIHOOD VALUES ARE FOUND*/
#define LSPNAME       100       /* # characters in sequence names */
#define MAXBRLENS 7.0 /* the maximum branch length of the species tree */
#define COLLAPSEBRLENS 1e-6 /* the branch length for collapsing branches */  
#define ERROR 1
#define NO_ERROR 0
#define YES 1
#define NO 0
#define NA -1
#define DEBUG 0

typedef struct sptable
	{
	char **taxaname;
	char **speciesname;
	int *taxanodenumber;
	}
	SPtable;

typedef struct node 
	{
	int father, nson, sons[NTAXA], namenumber;
	char taxaname[LSPNAME];
	double brlens, theta, support;
   	}
	Treenode;

typedef struct Tree
	{
   	int root;
	int ntaxa;
	int isrooted; 
	int numnodes;
	int outgroupnode;
	char outgroup[LSPNAME];
	int postorder[2*NTAXA];
	int preorder[2*NTAXA];	
   	Treenode nodes[2*NTAXA];
	}  
	Tree;


