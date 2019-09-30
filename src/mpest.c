/***********************************************************************
 *  mpest 2.0
 *
 *  copyright 2014-2019
 *
 *  Liang Liu
 *  Department of Statistics
 *  University of Georgia
 *
 *  lliu@uga.edu
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details (www.gnu.org).
 *
 * 
 ************************************************************************/

#include	"mpest.h"

int 		addNode (Tree *tree, int fromnode, int fromfather, int tonode);
int 		AddToPrintString (char *tempStr);
int 		Algorithm (int **triple, FILE *outfile, FILE *outputfile);
void 		copyTree(Tree *from, Tree *to);
int 		deleteNode (Tree *tree, int inode);
int 		genetreeTriples (Tree *tree, FILE *outfile);
int 		findNgenesandNtaxa (FILE *fTree);
int 		findNameNumber (Tree *tree);
int 		findTriple (int node1, int node2, int node3, long int *location);
void 		findOffsprings (int *offsprings, Tree *tree, int inode);
int 		findOutgroup (Tree *tree);
int 		logbinomialP (int n, int x, double p, double *logp);
int 		logLikelihood (int **triple, Tree *tree, double *loglike);
int 		maximizeaBrlen (int **triple, Tree *tree);
void 		MoveBrlens (Tree *tree);
int 		MoveNode (Tree *tree);
void		MrBayesPrint (char *format, ...);
void 		PrintHeader (void);
int 		printSptree (void);
int 		PrintState (int round, FILE *fout, int addend);
int 		PrintTree (Tree *tree, int inode, int showName, int showBrlens, int showTheta, int isRooted);
void 		printTriples (FILE *outfile);
void 		randomTree(Tree *tree);
void 		randomVector (int *array, int number);
int 		ReadaTree (FILE *fTree,Tree *tree);
void 		*SafeMalloc(size_t s);
int 		SaveSprintf(char **target, int *targetLen, char *fmt, ...);
int 		swapNodes (Tree *tree, int inode, int jnode);
int 		tripleDist (Tree *tree1, Tree *tree2);
int 		triples (Tree *tree, int ntrees);
int 		triplesInatree (int **triple, Tree *tree);
int			triplesInatree1	(int **triple, Tree *tree);
void 		WriteTreeToFile (Tree *tree, int inode, int showName, int showBrlens, int showTheta, int isRooted);
void		WritePhylipTreeToFile (Tree *tree, int inode, int showBrlens, int showTheta, int isRooted);
int			PrintPhylipTree(Tree *tree, int inode, int showBrlens, int showTheta, int isRooted);

Tree 		sptree;
Tree		gtree[NGENE];
int		    nGene = 0;
int		    **triplematrix; /*store the frequency of species tree triples among gene trees*/ 
char		spacer[10]="  ";
char		*printString;                /* string for printing to a file                */
size_t		printStringSize;             /* length of printString                        */
long int	seed = -1;
int 		usertree=0;		/* 1: use the usertree as the starting tree */
int		    taxanodenumber[NTAXA];
char		taxanames[NTAXA][LSPNAME];
int 		totaltaxa;
double		curLn;
int 		nruns;
int		    numupdatenodes;
int		    updatenodes[NTAXA];

int main (int argc, char *argv[])
{
	int i, j, k, n, ngene, distance; 
	FILE *fTree, *fout, *fin, *foutput, *foutput1, *foutput2;
	time_t t;
	struct tm *current;
	char genetreefile[LSPNAME], speciesname[NTAXA][LSPNAME], name[LSPNAME];
	char outfile[LSPNAME], outputfile[LSPNAME], outputfile1[LSPNAME], outputfile2[LSPNAME];

    PrintHeader();
	
    /*read the control file*/
	fin = (FILE*)gfopen(argv[1],"r");
	fscanf(fin,"%s%d%ld%d%d%d", genetreefile, &distance, &seed, &nruns, &ngene, &(sptree.ntaxa));
	sprintf(outfile, "%s_besttree.tre", genetreefile);
    sprintf(outputfile, "%s_output.tre", genetreefile);
	sprintf(outputfile1, "%s_genetree.triple", genetreefile);
	sprintf(outputfile2, "%s_genetree.dis", genetreefile);
	
    fout = (FILE*)gfopen(outfile,"w");
    foutput = (FILE*)gfopen(outputfile,"w");
	if(distance == 1){
		foutput1 = (FILE*)gfopen(outputfile1,"w");
		foutput2 = (FILE*)gfopen(outputfile2,"w");
	}
    fTree = (FILE*)gfopen(genetreefile,"r");
    
    /*set seed*/
	if(seed < 0){
		time(&t);
		current = localtime(&t);
		seed = 11*current->tm_hour + 111*current->tm_min + 1111*current->tm_sec + 372826; /*the random number generator requires that seed > 10000*/
		SetSeed(seed);
	}else{
        seed += 372826;
		SetSeed(seed);
    }
    
	/*the correspondence between species and individual sequences*/
    totaltaxa = 0;
	for(i=0; i<sptree.ntaxa; i++)
	{
        fscanf(fin,"%s%d", speciesname[i], &n);
		for(j=0; j<n; j++)
		{
			fscanf(fin,"%s", taxanames[totaltaxa]);	
			/*defines the species number that each taxon belongs to*/ 
			taxanodenumber[totaltaxa++] = i;
		}
	}	
     
    /*the lengths of the external branches are 1.0, because they are unestimable*/
    for(i=0; i<sptree.ntaxa; i++) sptree.nodes[i].brlens = 1.0;
    
    /****************************************************************************************
    read the user tree
    usertree=0: random starting tree
    usertree=1: use the usertree as the starting tree 
    usertree=2: only optimize branch lengths
    usertree=3: calculate likelihood score for a fixed tree
    usertree=4: 
    ****************************************************************************************/
	fscanf(fin,"%d", &usertree);
	if(usertree > 0){
		if(ReadaTree(fin, &sptree)== ERROR) {printf("Something is wrong with the user tree; It must be a rooted binary tree.\n");exit(-1);}
        
        /*we need to update taxanumber, because we change the order of the taxanames in the species tree*/
        for(i=0; i<totaltaxa; i++)
        {
            for(j=0; j<sptree.ntaxa; j++)
            {
                if(taxanodenumber[i] == j)
                {
                    for(k=0;k<sptree.ntaxa;k++)
                    {
                        if(strcmp(speciesname[j],sptree.nodes[k].taxaname)==0)  
                        {
                            taxanodenumber[i] = k;
                            break;
                        }
                    }
                    if(k == sptree.ntaxa)
                    {
                        printf("%s in the species-allele table is missing in the user-specified species tree\n\n",speciesname[j]);
                        exit(-1);
                    }
                    break;
                }
            }
        }
        
        /* read the taxa whose placements need to be updata*/	
        if(usertree == 4) 
	    {
		      fscanf(fin,"%d", &numupdatenodes);
		      for(j=0; j<numupdatenodes; j++)
		      {
 			        fscanf(fin,"%s",name);
			        for(k=0;k<sptree.ntaxa;k++)
			        {
				            if(strcmp(name,sptree.nodes[k].taxaname)==0)
                        	{
                            		updatenodes[j] = k;
                            		break;
                        	}
			         }
			         printf("%d \n",updatenodes[j]);
			         if(k == sptree.ntaxa)
			         {
				        printf("the taxa in the species tree do not have the update taxon %s\n",name);
				        exit(-1);
			         }
		      }
		      for(j=0;j<2*sptree.ntaxa-1;j++) if(sptree.nodes[j].brlens>0) sptree.nodes[j].brlens = 0.5;
	    }	
        
	}else{
        for(i=0; i<sptree.ntaxa; i++)
        {
            strcpy(sptree.nodes[i].taxaname, speciesname[i]);
        }
    }

	fclose (fin);

    /*find the number of gene trees in the gene tree file; mpest can use a subset of gene trees by setting a ngene (smaller than nGene) in the control file*/
	nGene = findNgenesandNtaxa(fTree);
	
	if(nGene > ngene) nGene = ngene;
	if(nGene < ngene){printf("Check on the number of gene trees in the control file. There are only %d gene trees in the tree file!\n",nGene);return(ERROR);}

	for(i=0; i<nGene; i++)
	{
		if(ReadaTree(fTree, &gtree[i])== ERROR) {printf("Something is wrong in the gene tree %d; It must be a rooted binary tree.\n",i+1);exit(-1);}
	}
	
	if(findNameNumber(gtree) == ERROR)
	{
		printf("Errors in findNameNumber function\n");
		exit(-1);
	}
	
	/*summarize triples in gene trees and the result is stored in triplematrix*/
	if(triples (gtree, nGene) == ERROR)
	{
		printf("Errors in the triples function\n");
		exit(-1);
	}
	
	if(distance == 1){
		/*calculate gene tree triples*/
		for(i=0; i<nGene; i++){
			fprintf(foutput1,"gene%d\t",i+1);
			genetreeTriples(&gtree[i], foutput1);
			fprintf(foutput1,"\n");
		}

		/*calculate triple distance*/
		for(i=0; i<nGene; i++){
			for(j=0; j<nGene; j++){
				if(j<i+1) fprintf(foutput2, "%d\t",0);
				else fprintf(foutput2, "%d\t", tripleDist(&gtree[i], &gtree[j]));
			}
			fprintf(foutput2,"\n");
		}
		return(NO_ERROR);
	}
	
    /*print the header of the output file*/
	if (PrintState (1, fout, 0) == ERROR)
        {
                MrBayesPrint ("%S   Errors in PrintState.\n", spacer);
                return (ERROR);
        }
    if (PrintState (1, foutput, 0) == ERROR)
    {
        MrBayesPrint ("%S   Errors in PrintState.\n", spacer);
        return (ERROR);
    }
    
    /*algorithm begins here*/
	for(i=0; i<nruns; i++)
	{
		if(Algorithm (triplematrix, fout, foutput) == ERROR)
		{
		printf("Errors in Algorithm\n");
		exit(-1);
		}
	}
	fprintf(fout, "end;\n\n");
    fprintf(foutput, "end;\n\n");

	/*free memory*/
	free(triplematrix[0]);
	free(triplematrix);
  	fclose(fTree);
	fclose(fout);
    fclose(foutput);
	if(distance == 1){
		fclose(foutput1);
		fclose(foutput2);
	}

  	return(NO_ERROR);
}

int logLikelihood (int **triple, Tree *tree, double *loglike)
{
	int i, j, k, w, son0, son1, father, cnode, ntriples=0;
	int offsprings0[NTAXA], offsprings1[NTAXA], offsprings2[NTAXA], taxa0[NTAXA], taxa1[NTAXA], taxa2[NTAXA], n0, n1, n2;
	long int location[2];
	double p, brlens, logp;
	int *array;

	if(DEBUG)
		array = (int*)calloc(sptree.ntaxa*(sptree.ntaxa-1)*(sptree.ntaxa-2)/6, sizeof (int));
 	
	*loglike = 0.0;
	for(i=tree->ntaxa; i<2*tree->ntaxa-1; i++)
	{
		if(i == tree->root)
			continue;
		else
		{
			brlens = tree->nodes[i].brlens;
			n0 = n1 = 0;
			for(j=0; j<tree->ntaxa; j++)
			{
				offsprings0[j] = 0;
				offsprings1[j] = 0;
			}

			son0 = tree->nodes[i].sons[0];
			son1 = tree->nodes[i].sons[1];
			findOffsprings (offsprings0, tree, son0);
			findOffsprings (offsprings1, tree, son1);

			for(j=0; j<tree->ntaxa; j++)
			{
				if (offsprings0[j] == 1)
					taxa0[n0++] = j;
				if (offsprings1[j] == 1)
					taxa1[n1++] = j;
			}

			cnode = i;
			father = tree->nodes[i].father;
			while (father != -1)
			{
				for(j=0; j<tree->ntaxa; j++)
					offsprings2[j] = 0;
						
				n2 = 0;
				p = 1-exp(-brlens)*2/3;

				if(tree->nodes[father].sons[0] == cnode)
					findOffsprings (offsprings2, tree, tree->nodes[father].sons[1]);
				else
					findOffsprings (offsprings2, tree, tree->nodes[father].sons[0]);

				for(j=0; j<tree->ntaxa; j++)
				{
					if (offsprings2[j] == 1)
						taxa2[n2++] = j;
				}

				for(j=0; j<n0; j++)
				{
					for(k=0; k<n1; k++)
						for(w=0; w<n2; w++)
						{
							findTriple (taxa0[j], taxa1[k], taxa2[w], location);
							if (DEBUG)
							{
								array[location[0]] = 1;
								/*printf("location %d %d %d %ld %ld %f %d %d\n",taxa0[j],taxa1[k],taxa2[w],location[0],location[1],brlens,triplematrix[location[0]][3],triplematrix[location[0]][location[1]]);*/
							}
							if(logbinomialP (triplematrix[location[0]][3],triplematrix[location[0]][location[1]],p,&logp) == ERROR)
							{
								printf("Errors in logbinomialP\n");
								return (ERROR);
							}
							/*printf("logp %f %d %d %d\n",logp, taxa0[j],taxa1[k],taxa2[w]);*/
							*loglike += logp;
						}
				}
				brlens += tree->nodes[father].brlens;
				cnode = father;
				father = tree->nodes[father].father;

				/*check the number of triples*/
				ntriples += n0 * n1 * n2;
			}					
		}
	}

	if(DEBUG)
	{
		for(i=0; i<sptree.ntaxa*(sptree.ntaxa-1)*(sptree.ntaxa-2)/6; i++)
			if(array[i] != 1)
			{
				printf("Errors in findTriple\n");
				return(ERROR);
			}
	}

	if(ntriples != (tree->ntaxa * (tree->ntaxa-1) * (tree->ntaxa-2) / 6))
		return (ERROR);

	return(NO_ERROR);
}

int logbinomialP (int n, int x, double p, double *logp)
{

	/**logp = 0.0;

	if(p>1.0 || p<0.0 || x > n)
	{
		printf("probability cannot be %f\n",p); 
		return(ERROR);
	}

	for(i=x+1; i<=n; i++)
		(*logp) += log(i);
	for(i=1; i<=(n-x); i++)
		(*logp) -= log(i);
	if(p == 1)
		(*logp) += x*log(0.999999) + (n-x)*log(1-0.999999);
	else
		(*logp) += x*log(p) + (n-x)*log(1-p);*/

	if(p == 1) p = 0.99999999;
	*logp = x*log(p) + (n-x)*log((1-p)/2);
	return (NO_ERROR);
}

int Algorithm (int **triple, FILE *outfile, FILE *outputfile)
{
	double oldloglike, diff;
	int i, round=1, sample = 1000, totalround = MAXROUND, stop = 0, noincreasing=0, accept = 0;
	Tree oldsptree;
	time_t starttime, endtime;
	clock_t start, end;
     	double cpu_time_used;

    /*starting tree*/
	if(usertree == 0) randomTree (&sptree);

	copyTree (&sptree, &oldsptree);
	for(i=0; i<sptree.ntaxa; i++) strcpy(oldsptree.nodes[i].taxaname, sptree.nodes[i].taxaname);
		
	if(logLikelihood (triple, &sptree, &curLn) == ERROR)
	{
		printf("Errors in loglikelihood\n");
		return (ERROR);
	}

    /*calculating the likelihood score for a fixed tree, which is given as the user tree*/
    if(usertree == 3) {printf("%lf\n",curLn);exit(-1);}

    /*record cpu time*/
	start = clock();
	starttime = time(0);

	while (round < totalround)
	{
		/*copy likelihood and tree*/
		copyTree (&sptree, &oldsptree);
		oldloglike = curLn;
		
		if(MoveNode (&sptree) == ERROR)
		{
			printf("Errors in MoveNode\n");
			return (ERROR);
		}
		if(logLikelihood (triple, &sptree, &curLn) == ERROR)
		{
			printf("Errors in logLikelihood\n");
			return (ERROR);
		}
		diff = curLn - oldloglike;

		if(diff < 0){
			noincreasing++;
			curLn = oldloglike;
			copyTree (&oldsptree, &sptree);
		}else{
			if(diff > 0)
				noincreasing = 0;			
			accept++;
		}
		
		round ++;
        if((round % sample) == 0 || stop == NUM_NOCHANGE){
            printf("\tround %d \t\tloglike %f \t....completed....\n",round, curLn);
            if (PrintState (round, outputfile, 0) == ERROR)
            {
                MrBayesPrint ("%S   Errors in PrintState.\n", spacer);
                return (ERROR);
            }
        }
		
		if(stop == NUM_NOCHANGE && round != 0)
		{
			if (PrintState (round, outfile, 0) == ERROR)
			{
                        MrBayesPrint ("%S   Errors in PrintState.\n", spacer);
                        return (ERROR);
			}
		}
		if (noincreasing > NUM_NOCHANGE || noincreasing == NUM_NOCHANGE) break;	
	}
	
	if (PrintState (round, outfile, 1) == ERROR)
        {
    		MrBayesPrint ("%S   Errors in PrintState.\n", spacer);
              	return (ERROR);
      	}

	/*record cpu time*/
	end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	endtime = time(0);
	if (difftime (endtime, starttime) > 2.0)
	{
		printf("\t\t\tAnalysis is completed in %.0f seconds (cpu time is %.0f seconds)\n", difftime(endtime, starttime), cpu_time_used);
	}
    else if (difftime (endtime, starttime) >= 1.0)
	{
		printf("\t\t\tAnalysis completed in 1 second\n");
	}
	else
	{
		printf("\t\t\tAnalysis completed in less than 1 second\n");
	}

	return (NO_ERROR);
}

void printTriples (FILE *outfile)
{
	int j, k;
	long int ntriples = sptree.ntaxa*(sptree.ntaxa-1)*(sptree.ntaxa-2)/6;

	fprintf(outfile, "  2. summary of the dataset\n");
	fprintf(outfile, "  \t\tntaxa\t\tnmissing\n");
	for(j=0;j<nGene;j++)
	{
		fprintf(outfile, "  gene%d\t\t%d\t\t%d\n", j, gtree[j].ntaxa, totaltaxa-gtree[j].ntaxa);
	}
	
	fprintf(outfile, "\n  3. triple matrix\n");
	fprintf(outfile, "  first\t\tsecond\t\tthird\t\tsum\n");
	for(j=0;j<ntriples;j++)
	{
		fprintf(outfile, "  ");
		for(k=0;k<4;k++)
			fprintf(outfile, "%d\t\t",triplematrix[j][k]);
		fprintf(outfile,"\n");
	}
}

int PrintState (int round, FILE *outfile, int addend)
{
	int i;
	char buffer[30];
  	struct timeval tv;
  	time_t curtime;

	/*print to file*/
	if (addend == 0)
	{
		if(round == 1)
    		{
				gettimeofday(&tv, NULL);
				curtime = tv.tv_sec;
				strftime(buffer,30,"%I:%M%p on %m-%d-%Y",localtime(&curtime));

				fprintf(outfile, "#Nexus\n[This mpest2.0 analysis was conducted at local time %s with seed = %ld. The branch length 1.0 of the external branches is arbitrary because mp-est does not estimate the lengths of the external branches. If the internal branch length is 7.0, it indicates that all gene trees support the species tree triple. In the output.tre file, the tree block contains the output trees (and their likelihood scores) generated from the algorithm, while the besttree.tre file contains the best tree. If multiple runs are specified, the besttree.tre file contains multiple best trees generated from multiple runs. Choose the one with the maximum likelihood score as the estimate of the species tree.]\n\n", buffer, seed);
				fprintf(outfile, "Begin trees;\n  translate\n");
				for (i=0; i<sptree.ntaxa-1; i++) {
					fprintf(outfile,"\t%d %s,\n", i+1, sptree.nodes[i].taxaname);
				}
				fprintf(outfile,"\t%d %s;\n", sptree.ntaxa, sptree.nodes[sptree.ntaxa-1].taxaname);
            }else{
                if (PrintTree(&sptree, sptree.root, 0, 1, 0, 1) == ERROR){
                    printf("Errors in printtree!\n");
                    return (ERROR);
                }
                fprintf(outfile, "  tree round%d [%2.6f] = %s", round, curLn, printString);
                free (printString);
            }
	}else{
		if (PrintTree(&sptree, sptree.root, 0, 1, 0, 1) == ERROR){
			printf("Errors in printtree!\n");
			return (ERROR);
		}
		fprintf(outfile, "  tree mpest [%2.6f] = %s", curLn, printString);
		free (printString);
	}
	return (NO_ERROR);
}

/*this function is for debugging purpose*/
int printSptree (void)
{
	int i;

    if (PrintTree(&sptree, sptree.root, 1, 1, 0, 1) == ERROR){
			printf("Errors in printtree!\n");
			return (ERROR);
	}
	printf("%s", printString);
	free (printString);

	for(i=0; i<2*sptree.ntaxa-1; i++)
		printf("node %d %s %d %d %d\n", i, sptree.nodes[i].taxaname, sptree.nodes[i].father,sptree.nodes[i].sons[0],sptree.nodes[i].sons[1]);

    return(NO_ERROR);
}

int MoveNode (Tree *tree)
{
	int inode, jnode, father, grandfather, nodes[5], nnodes=0, son;
	double rand = rndu();
	double p[2]={0.3,0.8};

	if(usertree == 2) {p[0] = 0.0; p[1]=0.0;}
    if(usertree == 4) {p[0] = 0.0; p[1] = 0.7;}
    
	if(rand < p[0])
	{
	   do{
		do{
			inode = (int)(rndu() * (2*tree->ntaxa-1));
		}while (inode == tree->root || inode == tree->nodes[tree->root].sons[0] || inode == tree->nodes[tree->root].sons[1]);

		father = tree->nodes[inode].father;
		grandfather = tree->nodes[father].father;

		if(tree->nodes[father].sons[0] == inode)
			son = tree->nodes[father].sons[1];
		else
			son = tree->nodes[father].sons[0];

		if(son >= tree->ntaxa)
		{
			nodes[nnodes++] = tree->nodes[son].sons[0];
			nodes[nnodes++] = tree->nodes[son].sons[1];
		}		
			
		if(tree->nodes[grandfather].sons[0] == father)
			nodes[nnodes++] = son = tree->nodes[grandfather].sons[1];
		else
			nodes[nnodes++] = son = tree->nodes[grandfather].sons[0];
	

		if(son >= tree->ntaxa)
		{
			nodes[nnodes++] = tree->nodes[son].sons[0];
			nodes[nnodes++] = tree->nodes[son].sons[1];
		}
	    }while(nnodes == 0);

		jnode = nodes[(int)(rndu() * nnodes)];
	
		swapNodes (tree, inode, jnode);
	}
	else if (rand < p[1])
	{
		do{
			inode = rndu() * tree->ntaxa;
            if(usertree == 4) inode = updatenodes[(int)(rndu()*numupdatenodes)];
		}while (tree->nodes[inode].father == tree->root);

        /*delete the branch of inode*/
	    if(deleteNode(tree, inode) == ERROR)
	    {
		 printf("Errors in deleteNode\n");
		 return (ERROR);
	    }

		do{
			jnode = rndu() * (2*tree->ntaxa-1);
		}while (jnode == tree->root || jnode == tree->nodes[inode].father || jnode == inode);

        /*attach inode to the branch of jnode*/
        if(addNode(tree, inode, tree->nodes[inode].father, jnode) == ERROR)
	    {
		  printf("Errors in addNode\n");
		  return (ERROR);
	    }
	}
	else
		MoveBrlens (tree);

	return (NO_ERROR);
}

void MoveBrlens (Tree *tree)
{
	int inode;
	double brlens, window = 0.1;

	do{
		inode = tree->ntaxa + (int)(rndu() * (tree->ntaxa-1));
	}while (inode == tree->root);

	brlens = tree->nodes[inode].brlens;
	brlens = brlens + (rndu()-0.5) * window;
	if(brlens < 0)
		brlens = 0.000001;
	if(brlens > MAXBRLENS)
		brlens = MAXBRLENS;
	tree->nodes[inode].brlens = brlens;
}

void copyTree(Tree *from, Tree *to)
{
   	int   i ;
 
  	/*copy the species tree nodes*/
  	to->root = from->root;
	to->ntaxa = from->ntaxa;
  	for(i=0; i<(2*from->ntaxa-1); i++)
  	{
    		to->nodes[i].nson = from->nodes[i].nson;
    		to->nodes[i].sons[0] = from->nodes[i].sons[0];
    		to->nodes[i].sons[1] = from->nodes[i].sons[1];
    		to->nodes[i].father = from->nodes[i].father;
    		to->nodes[i].brlens = from->nodes[i].brlens;
  	}
}

int findOutgroup (Tree *tree)
{
	int outgroup;

	if(tree->nodes[tree->root].sons[0] < tree->ntaxa)
		outgroup = tree->nodes[tree->root].sons[0];
	else
		outgroup = tree->nodes[tree->root].sons[1];
	return (outgroup);
}

int deleteNode (Tree *tree, int inode)
{
	int father, son, grandfather;
	
	if(inode == tree->root) return (ERROR);

	father = tree->nodes[inode].father;

	if(tree->nodes[father].sons[0] == inode) son = tree->nodes[father].sons[1];
	else son = tree->nodes[father].sons[0];

	if(father == tree->root)
	{
		tree->root = son;
	}
	else
	{
		grandfather = tree->nodes[father].father;
		if(tree->nodes[grandfather].sons[0] == father)
			tree->nodes[grandfather].sons[0] = son;
		else
			tree->nodes[grandfather].sons[1] = son;
		tree->nodes[son].father = grandfather;
	}
	return (NO_ERROR);
}
int addNode (Tree *tree, int fromnode, int fromfather, int tonode)
{
	int tofather;

    /*fromnode should not be the outgroup species*/
	if(fromfather == -1)
	{
		printf("Errors for fromfather\n");
		return (ERROR);
	}

    /*update the ancestral informatino of fromfather*/
	tofather = tree->nodes[tonode].father;
	if(tree->nodes[tofather].sons[0] == tonode)
		tree->nodes[tofather].sons[0] = fromfather;
	else
		tree->nodes[tofather].sons[1] = fromfather;
	tree->nodes[fromfather].father = tofather;

    /*update the ancestral information of tonode*/
	tree->nodes[tonode].father = fromfather;
	if(tree->nodes[fromfather].sons[0] == fromnode) tree->nodes[fromfather].sons[1] = tonode;
	else tree->nodes[fromfather].sons[0] = tonode;

	return (NO_ERROR);
}

int swapNodes (Tree *tree, int inode, int jnode)
{
	int ifather, jfather;
	ifather = tree->nodes[inode].father;
	jfather = tree->nodes[jnode].father;
	tree->nodes[inode].father = jfather;
	tree->nodes[jnode].father = ifather;
	if(tree->nodes[ifather].sons[0] == inode) 
		tree->nodes[ifather].sons[0] = jnode;
	else
		tree->nodes[ifather].sons[1] = jnode;
	if(tree->nodes[jfather].sons[0] == jnode)
		tree->nodes[jfather].sons[0] = inode;
	else
		tree->nodes[jfather].sons[1] = inode;
	return (NO_ERROR);
}
int maximizeaBrlen (int **triple, Tree *tree)
{
	return (NO_ERROR);
}
void randomVector (int *array, int number)
{
	int i, n, m;
	for(i=0; i<number; i++)
	{
		array[i] = i;
	}
	for(i=0; i<number; i++)
	{
		n = rndu() * (number - i);
		m = array[number-i-1];
		array[number-i-1] = array[n];
		array[n] = m;		
	}
}
void randomTree(Tree *tree)
{
	int i, outgroup, array[NTAXA];
	double brlens = 0.1;
	outgroup = gtree[0].nodes[findOutgroup (&gtree[0])].namenumber;
	randomVector (array, tree->ntaxa);
	tree->nodes[array[0]].nson = 0;
	tree->nodes[array[0]].sons[0] = -1;
	tree->nodes[array[0]].sons[1] = -1;
	tree->nodes[array[0]].father = tree->ntaxa;
	tree->nodes[tree->ntaxa].sons[0] = array[0];
	tree->nodes[array[0]].brlens = brlens;
	
	tree->nodes[array[1]].nson = 0;
	tree->nodes[array[1]].sons[0] = -1;
	tree->nodes[array[1]].sons[1] = -1;
	tree->nodes[array[1]].father = tree->ntaxa;
	tree->nodes[tree->ntaxa].sons[1] = array[1];
	tree->nodes[array[1]].brlens = brlens;
	
	tree->nodes[tree->ntaxa].father = tree->ntaxa + 1;
	tree->nodes[tree->ntaxa].brlens = brlens;
	tree->nodes[tree->ntaxa].nson = 2;
	
	tree->nodes[2*tree->ntaxa-2].nson = 2;
	tree->nodes[2*tree->ntaxa-2].sons[1] = 2*tree->ntaxa-3;
	tree->nodes[2*tree->ntaxa-2].father = -1;
	tree->nodes[2*tree->ntaxa-2].brlens = 0;
	
	for(i=2; i<tree->ntaxa; i++)
	{
		tree->nodes[array[i]].nson = 0;
		tree->nodes[array[i]].sons[0] = -1;
		tree->nodes[array[i]].sons[1] = -1;
		tree->nodes[array[i]].father = tree->ntaxa+i-1;
		tree->nodes[tree->ntaxa+i-1].sons[0] = array[i];
		tree->nodes[array[i]].brlens = brlens;
	}

	for(i=tree->ntaxa+1; i<2*tree->ntaxa-2; i++)
	{
		tree->nodes[i].nson = 2;
		tree->nodes[i].sons[1] = i-1;
		tree->nodes[i].father = i+1;
		tree->nodes[i].brlens = brlens;
	}
	tree->root = 2*tree->ntaxa-2;

	if(outgroup != array[tree->ntaxa - 1])
		swapNodes (tree, outgroup, array[tree->ntaxa-1]);
}


int triples (Tree *tree, int ntrees)
{
	int i, j, k;
	long int sum0=0, ntriples = sptree.ntaxa*(sptree.ntaxa-1)*(sptree.ntaxa-2)/6;

	for(i=0; i<ntrees; i++)
	{
		if(triplesInatree (triplematrix, &tree[i]) == ERROR)
		{
			printf("Errors in the tripleInatree function\n");
			return (ERROR);
		}
	}

	for(j=0;j<ntriples;j++)
	{
		for(k=0;k<3;k++)
			triplematrix[j][3] += triplematrix[j][k];
		sum0 += triplematrix[j][3];
	}
/*	for(j=0; j<nGene; j++)
		sum1 += (tree[j].ntaxa*(tree[j].ntaxa-1)*(tree[j].ntaxa-2)/6);

	if(sum0 != sum1)
	{
		printf("the number of triples %ld != %ld\n",sum0,sum1);
		return (ERROR);
	}
*/
	return (NO_ERROR);
}

void PrintHeader (void)
{

		MrBayesPrint ("\n\n\n\n%s            Maximum Pseudo-likelihood Estimation of Species Trees (MPEST_2.0)  \n\n",spacer);
		srand ((unsigned int)time(NULL));
			MrBayesPrint ("%s                                   by\n\n",spacer);
			MrBayesPrint ("%s                                Liang Liu\n\n",spacer);
			MrBayesPrint ("%s                        Department of Statistics\n",spacer);
			MrBayesPrint ("%s                          University of Georgia\n",spacer);
			MrBayesPrint ("%s                               lliu@uga.edu\n\n",spacer);
		    MrBayesPrint ("%s               Distributed under the GNU General Public License\n\n",spacer);	
}

int triplesInatree (int **triple, Tree *tree)
{
	int i, j, k, w, son0, son1;
	int offsprings0[NTAXA], offsprings1[NTAXA], offsprings2[NTAXA], n0, n1, n2;
	long int location[2];
 	
	
	for(i=tree->ntaxa; i<2*tree->ntaxa-1; i++)
	{
		if(i == tree->root)
			continue;
		else
		{
			n0 = n1 = n2 = 0;
			for(j=0; j<tree->ntaxa; j++)
			{
				offsprings0[j] = 0;
				offsprings1[j] = 0;
			}

			son0 = tree->nodes[i].sons[0];
			son1 = tree->nodes[i].sons[1];
			findOffsprings (offsprings0, tree, son0);
			findOffsprings (offsprings1, tree, son1);
			
			if(DEBUG)
			{
				printf("\n\n");
				for(j=0; j<tree->ntaxa; j++)
					printf("off %d %d %d %s %d %s\n", son0, son1, offsprings0[j],tree->nodes[offsprings0[j]].taxaname,offsprings1[j],tree->nodes[offsprings1[j]].taxaname);
			}
			

			for(j=0; j<tree->ntaxa; j++)
			{
				if(offsprings0[j] == 0 && offsprings1[j] == 0)
					offsprings2[n2++] = tree->nodes[j].namenumber;
				else if (offsprings0[j] == 1 && offsprings1[j] == 0)
					offsprings0[n0++] = tree->nodes[j].namenumber;
				else if (offsprings0[j] == 0 && offsprings1[j] == 1)
					offsprings1[n1++] = tree->nodes[j].namenumber;
				else
					return (ERROR);
			}

			for(j=0; j<n0; j++)
			{
				for(k=0; k<n1; k++)
				{
					if(offsprings0[j] == offsprings1[k])
						continue;

					for(w=0; w<n2; w++)
					{
						if(offsprings0[j] == offsprings2[w] || offsprings1[k] == offsprings2[w])
							continue;
						findTriple (offsprings0[j], offsprings1[k], offsprings2[w], location);
						triple[location[0]][location[1]]++;

						if(DEBUG)
						{
							printf("nodes %d %d %d location %ld %ld\n",offsprings0[j], offsprings1[k], offsprings2[w],location[0],location[1]);
						}

					}
				}
			}
					
		}
	}
	return(NO_ERROR);
}

int triplesInatree1 (int **triple, Tree *tree)
{
	int i, j, k, w, son0, son1, father, cnode;
	int offsprings0[NTAXA], offsprings1[NTAXA], offsprings2[NTAXA], n0, n1, n2;
	long int location[2];
	double brlens, minbrlens = COLLAPSEBRLENS;


	for(i=tree->ntaxa; i<2*tree->ntaxa-1; i++)
	{
		if(i == tree->root)
			continue;
		else
		{
			brlens = tree->nodes[i].brlens;
			n0 = n1 = 0;
			for(j=0; j<tree->ntaxa; j++)
			{
				offsprings0[j] = 0;
				offsprings1[j] = 0;
			}

			son0 = tree->nodes[i].sons[0];
			son1 = tree->nodes[i].sons[1];
			findOffsprings (offsprings0, tree, son0);
			findOffsprings (offsprings1, tree, son1);

			for(j=0; j<tree->ntaxa; j++)
			{
				if (offsprings0[j] == 1) offsprings0[n0++] = tree->nodes[j].namenumber;
				if (offsprings1[j] == 1) offsprings1[n1++] = tree->nodes[j].namenumber;
			}

			cnode = i;
			father = tree->nodes[i].father;
			
			while (father != -1)
			{	
				if(brlens > minbrlens){
					for(j=0; j<tree->ntaxa; j++) offsprings2[j] = 0;
							
					n2 = 0;

					if(tree->nodes[father].sons[0] == cnode)
						findOffsprings (offsprings2, tree, tree->nodes[father].sons[1]);
					else
						findOffsprings (offsprings2, tree, tree->nodes[father].sons[0]);

					for(j=0; j<tree->ntaxa; j++)
						if (offsprings2[j] == 1) 
							offsprings2[n2++] = tree->nodes[j].namenumber;

					for(j=0; j<n0; j++)
						for(k=0; k<n1; k++)
						{
							if(offsprings0[j] == offsprings1[k])
								continue;
							for(w=0; w<n2; w++)
							{
								if(offsprings0[j] == offsprings2[w] || offsprings1[k] == offsprings2[w])
									continue;
								findTriple (offsprings0[j], offsprings1[k], offsprings2[w], location);
								triple[location[0]][location[1]]++;
							}
						}
				}
				brlens += tree->nodes[father].brlens;
				cnode = father;
				father = tree->nodes[father].father;
			}					
		}
	}
	
	return(NO_ERROR);
}


void findOffsprings (int *offsprings, Tree *tree, int inode)
{
	int son0, son1;

	if(inode < tree->ntaxa)
		offsprings[inode] = 1;
	else
	{
		son0 = tree->nodes[inode].sons[0];
		son1 = tree->nodes[inode].sons[1];
		findOffsprings (offsprings, tree, son0);
		findOffsprings (offsprings, tree, son1);
	}
}

int findTriple (int n1, int n2, int n3, long int *location)
{
	int i, number=0, node1, node2, node3;

	if(n1 < n2 && n2 < n3)
	{
		node1 = n1;
		node2 = n2;
		node3 = n3;
		location[1] = 0;
	}
	else if(n2 < n1 && n1 < n3)
	{
		node1 = n2;
		node2 = n1;
		node3 = n3;
		location[1] = 0;
	}
	else if(n1 < n3 && n3 < n2)
	{
		node1 = n1;
		node2 = n3;
		node3 = n2;
		location[1] = 1;
	}
	else if(n2 < n3 && n3 < n1)
	{
		node1 = n2;
		node2 = n3;
		node3 = n1;
		location[1] = 1;
	}
	else if(n3 < n1 && n1 < n2)
	{
		node1 = n3;
		node2 = n1;
		node3 = n2;
		location[1] = 2;
	}
	else if(n3 < n2 && n2 < n1)
	{
		node1 = n3;
		node2 = n2;
		node3 = n1;
		location[1] = 2;
	}

	for(i=1; i<=node1; i++)
	{
		number += (sptree.ntaxa-i)*(sptree.ntaxa-i-1)/2;
	}
	for(i=node1+2; i<=node2; i++)
	{
		number += (sptree.ntaxa-i);
	}
	number += (node3 - node2 -1);
	location[0] = number;
	return(1);
}

int findNameNumber (Tree *tree)
{
	int i, j, k, stop;

	for(i=0; i<nGene; i++)
	{
		for(j=0; j<tree[i].ntaxa; j++)
		{
			stop = 0;
			for(k=0; k<totaltaxa; k++)
			{
				if(!strcmp(tree[i].nodes[j].taxaname, taxanames[k]))
				{
					stop = 1; 
					tree[i].nodes[j].namenumber = taxanodenumber[k];
					break;
				}
			}
			if(stop == 0)
			{
				printf("The taxon %s in gene tree %d is missing in the species-allele table in the control file\n",tree[i].nodes[j].taxaname, i);
				return (ERROR);
			}
		}
	}

	/*allocate memory for triplematrix*/
	triplematrix = (int**)calloc(sptree.ntaxa*(sptree.ntaxa-1)*(sptree.ntaxa-2)/6, sizeof(int*));
        triplematrix[0] = (int*)calloc(4*sptree.ntaxa*(sptree.ntaxa-1)*(sptree.ntaxa-2)/6, sizeof(int));
   	for(i = 0; i < (sptree.ntaxa*(sptree.ntaxa-1)*(sptree.ntaxa-2)/6); i++)
   		triplematrix[i] = triplematrix[0] + i*4;
  	if(!triplematrix)
	{
		printf(" allocating problem for triplematrix\n");
	   	return(ERROR);
	} 
	return(NO_ERROR);	
}

void MrBayesPrint (char *format, ...)
{
	va_list                 ptr;

	va_start (ptr, format);

			vprintf (format, ptr);
			fflush (stdout);

	va_end(ptr);

}

int findNgenesandNtaxa(FILE *fTree)
{
	char ch;
	int ncomma, ngene=0;

	ch = fgetc (fTree);
	while (ch != EOF) 
	{
      		ch = fgetc (fTree);
		ncomma = 0;
		while (ch != ';' && ch != EOF)
		{
			ch = fgetc (fTree);
			if(ch == ',')
				ncomma++;
		}
      		if(ch == ';')
		{
			gtree[ngene].ntaxa = ncomma + 1;
        		ngene++;
		}
	}
	rewind(fTree);
	return (ngene);
}

int ReadaTree (FILE *fTree,Tree *tree)
{
/* 
   Both names and numbers for species are accepted.  
   Species names are considered case-sensitive, with trailing blanks ignored.
*/
   	int cnode, cfather=-1, taxa = 0;  /* current node and father */
   	int inodeb=0;  /* node number that will have the next branch length */
   	int i, level=0, ch=' ';
   	char  skips[]="\"\'";
   	int nnode;   
	
	nnode = tree->ntaxa; 
   	FOR(i,2*(tree->ntaxa)-1) {
      		tree->nodes[i].father=-1;
      		tree->nodes[i].nson=0; 
			tree->nodes[i].sons[0] = -1;
			tree->nodes[i].sons[1] = -1;
   	}
	
   	while(ch != '(')
     	{
      		ch=fgetc(fTree);
     	}
   	ungetc(ch,fTree);

   	for (;;) {
      		ch = fgetc (fTree);
      		if (ch==EOF) return(-1);
      		else if (!isgraph(ch) || ch==skips[0] || ch==skips[1]) continue;
      		else if (ch=='(') {
         		level++;
         		cnode=nnode++;
   
         		if(nnode > 2*(tree->ntaxa)-1)
              		{
                  		printf("check tree: perhaps too many '('s");
                  		exit(-1);
              		}
         		if (cfather>=0) {
            			tree->nodes[cfather].sons[tree->nodes[cfather].nson++] = cnode;
            			tree->nodes[cnode].father=cfather;
         		}
         		else
            			tree->root=cnode;
         		cfather=cnode;
      		}
      		else if (ch==')') { level--;  inodeb=cfather; cfather=tree->nodes[cfather].father; }
      		else if (ch==':') fscanf(fTree,"%lf",&tree->nodes[inodeb].brlens);
      		else if (ch==',') ;
      		else if (ch==';' && level!=0) 
         	{
            		printf("; in treefile");
            		exit(-1);
         	}
      		else if (isdigit(ch))
      		{ 
         		ungetc(ch, fTree); 
         		fscanf(fTree,"%d",&inodeb); 
         		inodeb--;
         		tree->nodes[inodeb].father=cfather;
         		tree->nodes[cfather].sons[tree->nodes[cfather].nson++]=inodeb;
      		}
		else if (isalpha(ch))
		{		
			i = 0;
			while(ch != ':' && ch != ',' && ch != ')')
			{
				if(ch != ' ')
					tree->nodes[taxa].taxaname[i++] = ch;
				ch = fgetc(fTree);
			}
			ungetc(ch, fTree);
			tree->nodes[taxa].father = cfather;
			tree->nodes[cfather].sons[tree->nodes[cfather].nson++] = taxa;
			inodeb = taxa;
			taxa++;
		}
      		if (level<=0) break;
   	}
   
   	for ( ; ; ) {
      		while(isspace(ch=fgetc(fTree)) && ch!=';' );
      		if (ch==':')       fscanf(fTree, "%lf", &tree->nodes[tree->root].brlens);
      		else if (ch==';')  break;
      		else  { ungetc(ch,fTree);  break; }
   	}
	if(nnode != 2*(tree->ntaxa)-1) return(ERROR);
   
   	return (NO_ERROR);
}

/*PrinTree prints trees by taxa numbers; inode is the root of the tree*/
int PrintTree(Tree *tree, int inode, int showName, int showBrlens, int showTheta, int isRooted)
{

	char	*tempStr;
	int   tempStrSize;

	/* allocate the print string */
	printStringSize = 200;
	printString = (char *)SafeMalloc((size_t) (printStringSize * sizeof(char)));
	if (!printString)
		{
		MrBayesPrint ("%s   Problem allocating printString (%d)\n", spacer, printStringSize * sizeof(char));
		return (ERROR);
		}
	*printString = '\0';

	tempStrSize = 200;
	tempStr = (char *) SafeMalloc((size_t) (tempStrSize * sizeof(char)));
	if (!tempStr)
		{
		MrBayesPrint ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));
		return (ERROR);
		}

    /*print "(" and population size*/
	SaveSprintf (&tempStr, &tempStrSize,"(");
	AddToPrintString (tempStr);
					
	WriteTreeToFile (tree, tree->root, showName, showBrlens, showTheta, isRooted);

	if(showTheta == YES) 
		SaveSprintf (&tempStr, &tempStrSize,")[#%lf];\n",tree->nodes[tree->root].theta);
	else 
		SaveSprintf (&tempStr, &tempStrSize,");\n");
	AddToPrintString (tempStr);
	free (tempStr); 

	return (NO_ERROR);					
}

/*PrinPhylipTree prints trees by taxa names; inode is the root of the tree*/
int PrintPhylipTree(Tree *tree, int inode, int showBrlens, int showTheta, int isRooted)
{
	
	char			*tempStr;
	int                     tempStrSize;
	
	/* allocate the print string */
	printStringSize = 200;
	printString = (char *)SafeMalloc((size_t) (printStringSize * sizeof(char)));
	if (!printString)
	{
		MrBayesPrint ("%s   Problem allocating printString (%d)\n", spacer, printStringSize * sizeof(char));
		return (ERROR);
	}
	*printString = '\0';
	
	tempStrSize = 200;
	tempStr = (char *) SafeMalloc((size_t) (tempStrSize * sizeof(char)));
	if (!tempStr)
	{
		MrBayesPrint ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));
		return (ERROR);
	}
	
	SaveSprintf (&tempStr, &tempStrSize,"(");
	AddToPrintString (tempStr);
	
	WritePhylipTreeToFile (tree, tree->root, showBrlens, showTheta, isRooted);
	
	if(showTheta == YES) 
		SaveSprintf (&tempStr, &tempStrSize,")[#%lf];\n",tree->nodes[tree->root].theta);
	else 
		SaveSprintf (&tempStr, &tempStrSize,");\n");
	AddToPrintString (tempStr);
	free (tempStr); 
	
	return (NO_ERROR);					
	
}

void WriteTreeToFile (Tree *tree, int inode, int showName, int showBrlens, int showTheta, int isRooted)
{

		char			*tempStr;
		int      tempStrSize = 200;

		tempStr = (char *) SafeMalloc((size_t) (tempStrSize * sizeof(char)));
		if (!tempStr)
			MrBayesPrint ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));

		if (tree->nodes[inode].nson == 0)
		{
				if (showBrlens == YES){
                    if(showName == YES)SaveSprintf (&tempStr, &tempStrSize, "%s:1.0", tree->nodes[inode].taxaname);
    				else SaveSprintf (&tempStr, &tempStrSize, "%d:1.0", inode+1);
					AddToPrintString (tempStr);
				}else{
                    if(showName == YES)SaveSprintf (&tempStr, &tempStrSize, "%s", tree->nodes[inode].taxaname);
					else SaveSprintf (&tempStr, &tempStrSize, "%d", inode+1);
					AddToPrintString (tempStr);
				}
		}else{
				if (inode != tree->root)
				{
					SaveSprintf (&tempStr, &tempStrSize, "(");
					AddToPrintString (tempStr);
				}
				WriteTreeToFile (tree,tree->nodes[inode].sons[0],  showName, showBrlens, showTheta, isRooted);
				SaveSprintf (&tempStr, &tempStrSize, ",");
				AddToPrintString (tempStr);
				WriteTreeToFile (tree,tree->nodes[inode].sons[1], showName, showBrlens, showTheta, isRooted);	
				if (inode != tree->root)
				{
					if (tree->nodes[inode].father == tree->root && isRooted == NO)
					{
						if (showBrlens == YES)
						{
							SaveSprintf (&tempStr, &tempStrSize, ",%d:%lf", tree->nodes[inode].father + 1, tree->nodes[tree->nodes[inode].father].brlens);
							AddToPrintString (tempStr);
			
							if((tree->nodes[tree->nodes[inode].father].theta>0) && showTheta == YES) 
							{
								SaveSprintf (&tempStr, &tempStrSize, "[#%lf]", tree->nodes[tree->nodes[inode].father].theta);
								AddToPrintString (tempStr);
							}
						}
						else
						{
							SaveSprintf (&tempStr, &tempStrSize, ",%d", tree->nodes[inode].father + 1);
							AddToPrintString (tempStr);
						}
					}
				
					if (showBrlens == YES && isRooted == YES) /*tree->nodes[inode].father != tree->root)*/
					{
						SaveSprintf (&tempStr, &tempStrSize,"):%lf", tree->nodes[inode].brlens);
						AddToPrintString (tempStr);
						if((tree->nodes[inode].theta > 0) && showTheta == YES)
						{
							SaveSprintf (&tempStr, &tempStrSize, "[#%lf]", tree->nodes[inode].theta);
							AddToPrintString (tempStr);
						}
					}
					else
					{
						SaveSprintf (&tempStr, &tempStrSize, ")");
						AddToPrintString (tempStr);
					}					
		}
	}
	free (tempStr);	
}

void WritePhylipTreeToFile (Tree *tree, int inode, int showBrlens, int showTheta, int isRooted)
{
	
	char			*tempStr;
	int                      tempStrSize = 200;
	
	tempStr = (char *) SafeMalloc((size_t) (tempStrSize * sizeof(char)));
	if (!tempStr)
		MrBayesPrint ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));
	
	
	
	if (tree->nodes[inode].nson == 0)
	{
		/*if (showBrlens == YES)
		 {
		 SaveSprintf (&tempStr, &tempStrSize, "%d:%lf", inode+1, tree->nodes[inode].brlens);
		 AddToPrintString (tempStr);
		 if((tree->nodes[inode].theta>0) && showTheta == YES) 
		 {
		 SaveSprintf (&tempStr, &tempStrSize, "[#%lf]", tree->nodes[inode].theta);
		 AddToPrintString (tempStr);
		 }
		 }
		 else
		 {*/
		SaveSprintf (&tempStr, &tempStrSize, "%s:1.00", tree->nodes[inode].taxaname);
		AddToPrintString (tempStr);
		/*}*/
	}
	else
	{
		if (inode != tree->root)
		{
			SaveSprintf (&tempStr, &tempStrSize, "(");
			AddToPrintString (tempStr);
		}
		WritePhylipTreeToFile (tree,tree->nodes[inode].sons[0],  showBrlens, showTheta, isRooted);
		SaveSprintf (&tempStr, &tempStrSize, ",");
		AddToPrintString (tempStr);
		WritePhylipTreeToFile (tree,tree->nodes[inode].sons[1], showBrlens, showTheta, isRooted);	
		if (inode != tree->root)
		{
			if (tree->nodes[inode].father == tree->root && isRooted == NO)
			{
				if (showBrlens == YES)
				{
					SaveSprintf (&tempStr, &tempStrSize, ",%s:%lf", tree->nodes[tree->nodes[inode].father].taxaname, tree->nodes[tree->nodes[inode].father].brlens);
					AddToPrintString (tempStr);
					
					if((tree->nodes[tree->nodes[inode].father].theta>0) && showTheta == YES) 
					{
						SaveSprintf (&tempStr, &tempStrSize, "[#%lf]", tree->nodes[tree->nodes[inode].father].theta);
						AddToPrintString (tempStr);
					}
				}
				else
				{
					SaveSprintf (&tempStr, &tempStrSize, ",%s", tree->nodes[tree->nodes[inode].father].taxaname);
					AddToPrintString (tempStr);
				}
			}
			
			if (showBrlens == YES && isRooted == YES) /*tree->nodes[inode].father != tree->root)*/
			{
				SaveSprintf (&tempStr, &tempStrSize,"):%lf", tree->nodes[inode].brlens);
				AddToPrintString (tempStr);
				if((tree->nodes[inode].theta > 0) && showTheta == YES)
				{
					SaveSprintf (&tempStr, &tempStrSize, "[#%lf]", tree->nodes[inode].theta);
					AddToPrintString (tempStr);
				}
			}
			else
			{
				SaveSprintf (&tempStr, &tempStrSize, ")");
				AddToPrintString (tempStr);
			}					
		}
	}
	free (tempStr);
	
	
}

int genetreeTriples (Tree *tree, FILE *outfile)
{
	int **genetriplematrix;
	int i, j;
	long ntriples = sptree.ntaxa*(sptree.ntaxa-1)*(sptree.ntaxa-2)/6;

	/*allocate memory for genetriplematrix*/
	genetriplematrix = (int**)calloc(sptree.ntaxa*(sptree.ntaxa-1)*(sptree.ntaxa-2)/6, sizeof(int*));
        genetriplematrix[0] = (int*)calloc(4*sptree.ntaxa*(sptree.ntaxa-1)*(sptree.ntaxa-2)/6, sizeof(int));
   	for(i = 0; i < (sptree.ntaxa*(sptree.ntaxa-1)*(sptree.ntaxa-2)/6); i++)
   		genetriplematrix[i] = genetriplematrix[0] + i*4;
  	if(!genetriplematrix)
	{
		printf(" allocating problem for genetriplematrix1\n");
	   	return(ERROR);
	} 
	if(triplesInatree1 (genetriplematrix, tree) == ERROR)
	{
		printf("Errors in the tripleInatree function\n");
		return (ERROR);
	}	
	for(j=0;j<ntriples;j++)
	{
		if(genetriplematrix[j][0] == 1) fprintf(outfile,"1");
		else if(genetriplematrix[j][1] == 1) fprintf(outfile,"2");
		else if(genetriplematrix[j][2] == 1) fprintf(outfile,"3");
		else fprintf(outfile,"0");
		fprintf(outfile,"\t");	
	}
	
	free(genetriplematrix[0]);
	free(genetriplematrix);
	
	return(NO_ERROR);

}

int tripleDist (Tree *tree1, Tree *tree2)
{
	int **triplematrix1, **triplematrix2;
	int i, j, k;
	long ntriples = sptree.ntaxa*(sptree.ntaxa-1)*(sptree.ntaxa-2)/6;
	int dist = 0;

	/*allocate memory for triplematrix1*/
	triplematrix1 = (int**)calloc(sptree.ntaxa*(sptree.ntaxa-1)*(sptree.ntaxa-2)/6, sizeof(int*));
        triplematrix1[0] = (int*)calloc(4*sptree.ntaxa*(sptree.ntaxa-1)*(sptree.ntaxa-2)/6, sizeof(int));
   	for(i = 0; i < (sptree.ntaxa*(sptree.ntaxa-1)*(sptree.ntaxa-2)/6); i++)
   		triplematrix1[i] = triplematrix1[0] + i*4;
  	if(!triplematrix1)
	{
		printf(" allocating problem for triplematrix1\n");
	   	return(ERROR);
	} 

	/*allocate memory for triplematrix2*/
	triplematrix2 = (int**)calloc(sptree.ntaxa*(sptree.ntaxa-1)*(sptree.ntaxa-2)/6, sizeof(int*));
        triplematrix2[0] = (int*)calloc(4*sptree.ntaxa*(sptree.ntaxa-1)*(sptree.ntaxa-2)/6, sizeof(int));
   	for(i = 0; i < (sptree.ntaxa*(sptree.ntaxa-1)*(sptree.ntaxa-2)/6); i++)
   		triplematrix2[i] = triplematrix2[0] + i*4;
  	if(!triplematrix2)
	{
		printf(" allocating problem for triplematrix2\n");
	   	return(ERROR);
	} 

	if(triplesInatree1 (triplematrix2, tree2) == ERROR)
	{
		printf("Errors in the tripleInatree function\n");
		return (ERROR);
	}

	if(triplesInatree1 (triplematrix1, tree1) == ERROR)
	{
		printf("Errors in the tripleInatree function\n");
		return (ERROR);
	}

	for(j=0;j<ntriples;j++)	
		for(k=0; k<3; k++)
	 		dist += (triplematrix1[j][k]-triplematrix2[j][k]) * (triplematrix1[j][k]-triplematrix2[j][k]);
		
	free(triplematrix1[0]);
	free(triplematrix1);
	free(triplematrix2[0]);
	free(triplematrix2);

	return(dist/2);
}


int AddToPrintString (char *tempStr)

{

        size_t                  len1, len2;

        len1 = (int) strlen(printString);
        len2 = (int) strlen(tempStr);
        if (len1 + len2 + 5 > printStringSize)
                {
                printStringSize += len1 + len2 - printStringSize + 200;
                printString = realloc((void *)printString, printStringSize * sizeof(char));
                if (!printString)
                        {
                        MrBayesPrint ("%s   Problem reallocating printString (%d)\n", spacer, printStringSize * sizeof(char));
                        goto errorExit;
                        }
                }
        strcat(printString, tempStr);
#       if 0
        printf ("printString(%d) -> \"%s\"\n", printStringSize, printString);
#       endif
        return (NO_ERROR);

        errorExit:
                return (ERROR);

}



#define TARGETLENDELTA (100)

int SaveSprintf(char **target, int *targetLen, char *fmt, ...) {
  va_list    argp;
  int        len,retval;

  va_start(argp, fmt);
#ifdef VISUAL
  len = _vsnprintf(NULL, 0, fmt, argp);
#else
  len = vsnprintf(NULL, 0, fmt, argp);
#endif

  va_end(argp);

  if(len>*targetLen)
        {
/*        fprintf(stderr, "* increasing buffer from %d to %d bytes\n", *targetLen, len+TARGETLENDELTA); */
        *targetLen = len+TARGETLENDELTA; /* make it a little bigger */
            *target = realloc(*target, *targetLen);
        }

  va_start(argp,fmt);
  retval=vsprintf(*target, fmt, argp);
  va_end(argp);

/*   fprintf(stderr, "* savesprintf /%s/\n",*target); */
  return retval;
}

void *SafeMalloc(size_t s) {
        void *ptr = malloc(s);
        if(ptr==NULL)
                return NULL;
        return memset(ptr,0,s);
}



























