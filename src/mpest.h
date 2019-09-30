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

#define NTAXA         400      /* max # of species */
#define NGENE         20000      /* max # of loci */
#define MAXROUND	10000000		/* MAX # OF ROUNDS*/
#define NUM_NOCHANGE	20000		/* # OF ROUNDS THAT NO BIGGER LIKELIHOOD VALUES ARE FOUND*/
#define LSPNAME       60       /* # characters in sequence names */
#define MAXBRLENS 7.0 /* the maximum branch length of the species tree */
#define COLLAPSEBRLENS 1e-6 /* the branch length for collapsing branches */  
#define ERROR 1
#define NO_ERROR 0
#define YES 1
#define NO 0
#define NA -1
#define DEBUG 0
#define FPN(file) fputc('\n', file)
#define FOR(i,n) for(i=0; i<n; i++)
#define PointGamma(prob,alpha,beta) PointChi2(prob,2.0*(alpha))/(2.0*(beta))
#define CDFGamma(x,alpha,beta) IncompleteGamma((beta)*(x),alpha,LnGamma(alpha))
typedef struct node 
	{
	int father, nson, sons[2], namenumber;
	char taxaname[LSPNAME];
	double brlens,theta;
   	}
	Treenode;
typedef struct Tree
	{
   	int root;
	int ntaxa; 
   	Treenode nodes[2*NTAXA];
	}  
	Tree;
/* tool functions*/
FILE *gfopen(char *filename, char *mode);
void SetSeed (unsigned int seed);
double LnGamma (double x);
double rndu (void);
