/***********************************************************************
 *  mpest 3.0
 *
 *  copyright 2014-2024
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details (www.gnu.org).
 *
 * 
 ************************************************************************/

#include	"mpest.h"
#include 	"neighbour_joining.h"
#include 	"tool.h"

int 		AddNode (Tree *tree, int fromnode, int fromfather, int tonode);
int 		AllocateMemory (int nspecies);
int 		Algorithm (int **triple, FILE *besttreefile, FILE *outputfile, int runindex);
void 		CopyTree (Tree *from, Tree *to);
int 		DeleteNode (Tree *tree, int inode);
int 		FindNodeAncestors (int node, int *ancestors, Tree *tree);
int 		FindTriple (int node1, int node2, int node3, int *location);
int 		FindQuartet (int node1, int node2, int node3, int node4, int *location);
void 		FindOffsprings (int *offsprings, Tree *tree, int inode);
int 		FindOutgroup (Tree *tree);
void 		FreeMemory (void);
int 		GenetreePartitions (char *treefile, int ntrees);
int 		GenetreeNodeDistance (char *treefile, double *nodedistance);
int 		InitialSpeciesTreePara (void);
int 		LogBinomialP (int n, int x, double p, double *logp);
int 		Loglikelihood (int **triple, Tree *tree, double *loglike);
void 		MoveBrlens (Tree *tree);
int 		MoveNode (Tree *tree);
int 		NJst (double *nodedistance, char **speciesname, uint32_t species_count, int outgroupSpecies);
int 		NodeDistance (int node1, int node2, Tree *tree);
int 		NodeDistanceLocation (int sp1, int sp2, int *location);
void 		PrintHeader (void);
void 		PrintHelp (void);
int			PrintPhylipTree (Tree *tree, int inode, int showBrlens, int showTheta, int isRooted);
int 		printSptree (void);
int 		PrintState (int round, FILE *fout, int addend);
int 		PrintTree (Tree *tree, int inode, int showName, int showBrlens, int showTheta, int showSupport, int isRooted);
int 		QuartetDistance (char *quartetfile, int ntrees);
int 		QuartetsFreqInatree (int **quartet, Tree *tree);
int 		QuartetsList (char *treefile, int ntrees);
void 		RandomTree (Tree *tree);
void 		RandomVector (int *array, int number);
int 		ReadaTree (FILE *fTree, Tree *tree);
int 		ReadaTreeString (char *treestring, Tree *tree, bool taxaorder);
int 		ReadSPtableFromInputFile (char *inputfile);
int 		ReadParametersFromInputFile (char *inputfile);
int 		SptreePartitionSupport (char *treefile);
int 		StartSpeciesTree (int userdefinedtree);
int 		SwapNodes (Tree *tree, int inode, int jnode);
int 		TreePartitions (int **partitions, Tree *tree);
int			TreeTraversalPostorder (int node, Tree *tree);
int			TreeTraversalPreorder (int node, Tree *tree);
int 		TreeBranchCollapse (char *treefile);
int 		TripleDistance (char *treefile, int ntrees);
int 		TriplesFreq (char *treefile, int ngene);
int 		TriplesFreqInatree (int **triple, Tree *tree);
int			TriplesFreqInatreeCollapse (int **triple, Tree *tree); /*we do NOT need this function*/
int 		TriplesList (char *treefile, int ntrees);
int 		UpdateSPtable ();
void 		WriteTreeToFile (Tree *tree, int inode, int showName, int showBrlens, int showTheta, int showSupport, int isRooted);
void		WritePhylipTreeToFile (Tree *tree, int inode, int showBrlens, int showTheta, int isRooted);

Tree 		sptree;
int		    **triplematrix; /*store the frequency of gene tree triples*/ 
char		spacer[10]="  ";
long int	seed = -1;
int 		usertree=0;		/* 1: use the usertree as the starting tree */
char 		genetreefile[LSPNAME];
char		usertreefile[LSPNAME];
int 		totaltaxa;
double		curLn;
float 		collapseBrlens;
int 		nruns=1;
int		    numupdatenodes;
int			numGenes;
int			genetreetype; /*0: rooted 1:unrooted*/
int		    updatenodes[NTAXA];
int 		postorderindex = 0;
int			calculatedistance = 0;
SPtable		sptable;

int main (int argc, char *argv[]){
	int i, opt; 
	char besttreefile[150], outputtreefile[150];
	time_t t;
	struct tm *current;
	FILE *fbesttree, *foutputtree;

	/*parse command line*/
	while((opt = getopt(argc, argv, ":i:s:u:n:c:hTQPLBN")) != -1) { 
        switch(opt) { 
            case 'i': 
				strcpy(genetreefile, optarg); 
                break; 
            case 's': 
				if(atoi(optarg)){
					seed = atoi(optarg);
				}            
                break; 
			case 'u':
				strcpy(usertreefile, optarg);
				usertree = 1;
				break;
			case 'n':
				nruns = atoi(optarg);
				break;
			case 'c':
				calculatedistance = 4;
				collapseBrlens = atof(optarg);
				break;
			case 'h':
				PrintHelp();
				return NO_ERROR;
				break;
			case 'T':
				calculatedistance = 1;
				break;
			case 'Q':
				calculatedistance = 2;
				break;
			case 'P':
				calculatedistance = 3;
				break;
			case 'N':
				calculatedistance = 5;
				break;
			case 'L':
				usertree = 2;
				break;
			case 'B':
				usertree = 3;
				break;
            case ':': 
				PrintHelp();
                printf("Command line error: option -%c needs a value\n", optopt);
				return NO_ERROR; 
                break; 
            case '?': 
				PrintHelp();
                printf("Command line error: unknown option: -%c\n", optopt);
				return NO_ERROR;
                break; 
        } 
    } 

	if (genetreefile[0] == '\0'){
		PrintHelp();
		return ERROR;
	}

	if(usertree>0 && strlen(usertreefile) == 0){
		PrintHelp();
		if(usertree == 2){
			printf("Command line error: please provide a user tree for option -L\n");
		}
		if(usertree == 3){
			printf("Command line error: please provide a user tree for option -B");
		}
		return NO_ERROR;
	}

	/*read parameters from input file*/
	if(ReadParametersFromInputFile(genetreefile) == ERROR){
		printf("Errors in ReadParametersFromInputFile\n");
		return ERROR;
	}

	if(genetreetype == 1 && calculatedistance == 0){
		printf("\nThe input gene trees for the MPEST analysis must be rooted trees\nNote that input gene trees do not have to be binary. MPEST can take polytomy (rooted) gene trees.");
		return ERROR;
	}
	if(genetreetype == 1 && calculatedistance == 1){
		printf("\nTo calculate Triple Distance, input gene trees must be rooted trees.\nNote that polytomy (rooted) gene trees are allowed.");
		return ERROR;
	}
	
	/*allocate memory*/
	AllocateMemory (sptree.ntaxa);

	/*read sptable from input file*/
	if(ReadSPtableFromInputFile(genetreefile) == ERROR){
		printf("Errors in ReadSPtableFromInputFile\n");
		return ERROR;
	}

	/*initialize species tree*/
	if(InitialSpeciesTreePara() == ERROR){
		printf("Errors in InitialSpeciesTreePara\n");
		return ERROR;
	}

	/*set seed*/
	if(seed < 0){
		time(&t);
		current = localtime(&t);
		seed = 11*current->tm_hour + 111*current->tm_min + 1111*current->tm_sec;
		SetSeed(seed);
	}else{
		SetSeed(seed);
    }

	/************************************************************************
	*	 Algorithm for triple, quartet, partition, collapse, njst          *
	************************************************************************/
	if(calculatedistance == 1){
		if(TriplesList(genetreefile, numGenes) == ERROR){
			printf("Errors in TriplesList\n");
			return ERROR;
		}
		if(TripleDistance (genetreefile, numGenes) == ERROR){
			printf("Errors in TripleDistance\n");
			return ERROR;
		}
		printf("\nGene tree triples are saved in file %s_triple.list", genetreefile);
		printf("\nTriple distances among gene trees are saved in file %s_triple.dis\n", genetreefile);
		return NO_ERROR;
	}
	if(calculatedistance == 2){
		if(QuartetsList(genetreefile, numGenes) == ERROR){
			printf("Errors in QuartetsList\n");
			return ERROR;
		}
		if(QuartetDistance (genetreefile, numGenes) == ERROR){
			printf("Errors in QuartetDistance\n");
			return ERROR;
		}
		printf("\nGene tree quartets are saved in file %s_triple.list", genetreefile);
		printf("\nQuartet distances among gene trees are saved in file %s_triple.dis\n", genetreefile);
		return NO_ERROR;
	} 
	if(calculatedistance == 3){
		if(GenetreePartitions(genetreefile, numGenes) == ERROR){
			printf("Errors in GenetreePartitions\n");
			return ERROR;
		}
		printf("\nPartitions of gene trees are saved in file %s_partition\n", genetreefile);
		return NO_ERROR;
	} 
	if(calculatedistance == 4){
		if(TreeBranchCollapse(genetreefile) == ERROR){
			printf("Errors in TreeBranchCollapse\n");
			return ERROR;
		}
		printf("\nPolytomy gene trees are saved in file %s_collapse.tre\n", genetreefile);
		return NO_ERROR;
	}
	if(calculatedistance == 5){
		/*build NJst tree and use it as the start tree in the algorithm*/
		double *nodedistance;
		nodedistance = (double*) malloc (sptree.ntaxa * (sptree.ntaxa-1) / 2 * sizeof(double*));

		/*calculate nodedistance*/
		if(GenetreeNodeDistance (genetreefile, nodedistance) == ERROR){
			printf("Errors in GenetreeNodeDistance\n");
			return ERROR;
		}	

		/*build NJst tree*/
		NJst(nodedistance, sptable.speciesname, sptree.ntaxa, sptree.outgroupnode);

		/*set species tree*/
		ReadaTreeString (printString, &sptree, YES);
		sprintf(besttreefile, "%s_njst.tre", genetreefile);
		fbesttree = fopen(besttreefile,"w");
		fprintf(fbesttree,"begin trees;\n   Tree njst = %s\nEnd;", printString);
		printf("\nNJst tree is saved in file %s\n", besttreefile);
		
        free(printString);
		free(nodedistance);
		fclose(fbesttree);
		return NO_ERROR;
	}
	
	/************************************************************************
	*	Algorithm for MP-EST                                                *
	************************************************************************/
	PrintHeader();
	printf("        Build NJst tree (initial tree in algorithm)     ....");
	if(StartSpeciesTree (usertree) == ERROR){
		printf("Errors in StartSpeciesTree\n");
		return ERROR;
	};
	printf("completed....\n");
        
	/*summarize triples in gene trees*/
	printf("        Calculate triple frequencies of %d gene trees ....", numGenes);
	if(TriplesFreq (genetreefile, numGenes) == ERROR){
		printf("Errors in TriplesFreq\n");
		return ERROR;
	}
	printf("completed....\n\n        MPEST algorithm starts here....\n\n");

    /*algorithm begins here*/
	sprintf(besttreefile, "%s_besttree.tre", genetreefile);
	sprintf(outputtreefile, "%s_output.tre", genetreefile);
	fbesttree = fopen(besttreefile,"w");
	foutputtree = fopen(outputtreefile,"w");

	for(i=0; i<nruns; i++){
		if(Algorithm (triplematrix, fbesttree, foutputtree, i) == ERROR){
			printf("Errors in Algorithm\n");
			return ERROR;
		}
	}
	fprintf(fbesttree, "end;\n\n");
    fprintf(foutputtree, "end;\n\n");
	fclose(fbesttree);
	fclose(foutputtree);

	/*free memory*/
	FreeMemory();

  	return NO_ERROR;
}

void PrintHelp (void){
	char message[]= "\n\nMaximum Pseudo-likelihood Estimation of Species Trees - MPEST\nVersion 3.0\n(c) Copyright, 2014-2024 Liang Liu\nDepartment of Statistics, University of Georgia\n\nUsage: mpest [-i inputfile] [-n #] [-s #] [-u NAME] [-h] [-B|-C|-L|-N|-P|-Q|-T|]\n  -i: NAME = name of input file in nexus format. Input gene trees must be rooted; polytomy trees are allowed.\n  -n: # = number of runs [default = 1].\n  -s: # = seed for random number generator [default = system generated].\n  -u: NAME = name of user tree file.\n  -c: # = convert short branches (<#) to polytomies in gene trees.\n  -h: help message [default = random].\n  -B: optimize branch lengths of a fixed species tree provided through usertree option -u.\n  -L: calculate loglikelihood of a species tree provided through usertree option -u.\n  -N: build NJst tree.\n  -P: calculate partitions for gene trees.\n  -Q: calculate pairwise quartet distances among gene trees. Polytomy (unrooted or rooted) gene trees are allowed.\n  -T: calculate pairwise triple distances of gene trees. Polytomy rooted gene trees are allowed.\n\n";
  	printf("%s", message);
}

void FreeMemory (void){
	/*free memory*/
	free(triplematrix[0]);
	free(triplematrix);

	/*free sptable*/
	free(sptable.taxanodenumber);

	/*free sptable.speciesname*/	
    free(sptable.speciesname[0]);
	free(sptable.speciesname);

	/*free sptable.taxaname*/
	free(sptable.taxaname[0]);
	free(sptable.taxaname);
}

int StartSpeciesTree (int userdefinedtree){
	int j;
	FILE *fin;

	/****************************************************************************************
    read the user tree
    usertree=0: random starting tree
    usertree=1: use the usertree as the starting tree 
    usertree=2: only optimize branch lengths
    usertree=3: calculate likelihood score for a fixed tree
    usertree=4: 
    ****************************************************************************************/
	
	if(userdefinedtree == 0){
		/*build NJst tree and use it as the start tree in the algorithm*/
		double *nodedistance;
		nodedistance = (double*) malloc (sptree.ntaxa * (sptree.ntaxa-1) / 2 * sizeof(double*));

		/*calculate nodedistance*/
		if(GenetreeNodeDistance (genetreefile, nodedistance) == ERROR){
			printf("Errors in GenetreeNodeDistance\n");
			return ERROR;
		}	

		/*build NJst tree*/
		NJst(nodedistance, sptable.speciesname, sptree.ntaxa, sptree.outgroupnode);

		/*set species tree*/
		ReadaTreeString (printString, &sptree, YES);
        free(printString);

		/*if taxaorder=YES when setting sptree, no need to update SPtable*/
		//UpdateSPtable();

		free(nodedistance);
	}else{
		/*read user species tree*/
		fin = fopen(usertreefile, "r");
		if(!fin){
			printf("Error: cannot open file %s\n", usertreefile);
			return ERROR;
		}

		if(ReadaTree(fin, &sptree) == ERROR){
			printf("Errors in the user tree; It must be a rooted binary tree.\n");
			return ERROR;
		}
		fclose(fin);

        /*update SPtable*/
        UpdateSPtable();
		
        /*read taxa whose placements need to be updated*/	
       /* if(userdefinedtree == 4) {
			if(fscanf(fin,"%d", &numupdatenodes) != 1){
				printf("There are something wrong in the control file\n");
				return ERROR;
			}	
			for(j=0; j<numupdatenodes; j++){
				if(fscanf(fin,"%s",name) != 1){
					printf("There are something wrong in the control file\n");
					return ERROR;
				}	
				for(k=0;k<sptree.ntaxa;k++){
					if(strcmp(name,sptree.nodes[k].taxaname)==0){
						updatenodes[j] = k;
						break;
					}
				}
				printf("%d \n",updatenodes[j]);
				if(k == sptree.ntaxa){
					printf("the taxa in the species tree do not have the update taxon %s\n",name);
					exit(-1);
				}
			}
			for(j=0;j<2*sptree.ntaxa-1;j++) if(sptree.nodes[j].brlens>0) sptree.nodes[j].brlens = 0.5;
		}	*/    
	}

	/*terminal branch length = 1 if one allele per species*/
	if(totaltaxa == sptree.ntaxa){
		for(j=0; j<sptree.ntaxa; j++){
			sptree.nodes[j].brlens = 1.0;
		}
	}

	return NO_ERROR;
}

int AllocateMemory (int nspecies){
	int ntriples, i;

	/*allocate triplematrix*/
	ntriples = nspecies*(nspecies-1)*(nspecies-2)/6;
	triplematrix = (int**)calloc(ntriples, sizeof(int*));
    triplematrix[0] = (int*)calloc(4*ntriples, sizeof(int));
   	for(i=0; i<ntriples; i++)
   		triplematrix[i] = triplematrix[0] + i*4;
  	if(!triplematrix){
		printf(" allocating problem for triplematrix\n");
	   	return ERROR;
	} 

	/*allocate sptable*/
	sptable.taxanodenumber = (int*) malloc (totaltaxa * sizeof(int*));

	/*allocate sptable.speciesname*/
	sptable.speciesname = (char**)calloc(nspecies, sizeof(char*));
    sptable.speciesname[0] = (char*)calloc(nspecies*LSPNAME, sizeof(int));
   	for(i=0; i<nspecies; i++){
		sptable.speciesname[i] = sptable.speciesname[0] + i*LSPNAME;
	}
  	if(!sptable.speciesname){
		printf(" allocating problem for sptable.speciesname\n");
	   	return ERROR;
	} 

	/*allocate sptable.taxaname*/
	sptable.taxaname = (char**)calloc(totaltaxa, sizeof(char*));
    sptable.taxaname[0] = (char*)calloc(totaltaxa*LSPNAME, sizeof(int));
   	for(i=0; i<totaltaxa; i++){
		sptable.taxaname[i] = sptable.taxaname[0] + i*LSPNAME;
	}  		
  	if(!sptable.taxaname){
		printf(" allocating problem for sptable.taxaname\n");
	   	return ERROR;
	} 

	return NO_ERROR;
}

int ReadParametersFromInputFile (char *inputfile){
	int c, index, i, j, n;
	char buffer[LSPNAME], x[LSPNAME];	
	FILE *fin = fopen(inputfile, "r");

	if(!fin){
		PrintHelp();
		printf("Error: cannot open file %s\n", inputfile);
		return ERROR;
	}

	while(!feof(fin)){
		if(fscanf(fin, "%s", buffer) != 1){
			printf("Error: cannot read gene tree file 1 %s\n", inputfile);
			return ERROR;
		}
		
		for(char *p=buffer; *p; p++){
			*p=tolower(*p);
		} 

		/*sptree.ntaxa and numGenes*/
		if(!strcmp(buffer,"dimension")){
			c=fgetc(fin);
			while(c != 'n'){
				c=fgetc(fin);
			}
			c=tolower(fgetc(fin));
			if(c=='t'){
				while(c != '='){
					c=fgetc(fin);
				}
				if(fscanf(fin,"%d", &sptree.ntaxa) != 1){
					printf("Error: cannot read gene tree file 2 %s\n", inputfile);
					return ERROR;;
				}
				c=fgetc(fin);
				while(c != '='){
					c=fgetc(fin);
				}
				
				if(fscanf(fin,"%d", &numGenes) != 1){
					printf("Error: cannot read gene tree file 3 %s\n", inputfile);
					return ERROR;
				}
			}else{
				while(c != '='){
					c=fgetc(fin);
				}
				if(fscanf(fin,"%d", &numGenes) != 1){
					printf("Error: cannot read gene tree file %s\n", inputfile);
					return ERROR;;
				}
				c=fgetc(fin);
				while(c != '='){
					c=fgetc(fin);
				}
				if(fscanf(fin,"%d", &sptree.ntaxa) != 1){
					printf("Error: cannot read gene tree file 4 %s\n", inputfile);
					return ERROR;;
				}
			}
		}
		
		/*sptree.outgroup and treetype*/
		if(!strcmp(buffer,"format")){
			c=fgetc(fin);
			while(!isalpha(c)){
				c=fgetc(fin);
			}
			if(tolower(c)=='t'){
				while(c != '='){
					c=fgetc(fin);
				}
				if(fscanf(fin,"%s", x) != 1){
					printf("Error: cannot read gene tree file 5 %s\n", inputfile);
					return ERROR;;
				}
				if(tolower(x[0]) == 'r'){
					genetreetype = 0;
				}else{
					genetreetype = 1;
				}
				c = fgetc(fin);
				while(c != '='){
					c = fgetc(fin);
				}
				while(!isalpha(c)){
					c = fgetc(fin);
				}
				sptree.outgroup[0]=c;
				index = 1;
				while(c != ' ' && c != ';' && c != '\t'){
					c = fgetc(fin);
					sptree.outgroup[index++] = c;
				}
				sptree.outgroup[index-1] = '\0';
			}else{
				while(c != '='){
					c=fgetc(fin);
				}
				if(fscanf(fin,"%s", sptree.outgroup) != 1){
					printf("Error: cannot read gene tree file 6 %s\n", inputfile);
					return ERROR;
				}
				c = fgetc(fin);
				while(c != '='){
					c = fgetc(fin);
				}
				if(fscanf(fin,"%s", x) != 1){
					printf("Error: cannot read gene tree file 7 %s\n", inputfile);
					return ERROR;
				}
				if(tolower(x[0]) == 'r'){
					genetreetype = 0;
				}else{
					genetreetype = 1;
				}
			}			
		}

		/*find totaltaxa*/
		if(!strcmp(buffer,"matrix")){
			totaltaxa = 0;
			for(i=0; i<sptree.ntaxa; i++){
				if(fscanf(fin,"%s %d", x, &n) != 2){
					printf("# of species in SPtable != ntaxa in file %s\n", inputfile);
					return ERROR;
				}
				totaltaxa += n;
				for(j=0; j<n; j++){
					if(fscanf(fin,"%s", x) != 1){
						printf("Errors in SPtable in file %s\n", inputfile);
						return ERROR;
					}	
				}
			}			
		}
		
		if(!strcmp(buffer,"end;")){
			break;
		}
	}

	while(!feof(fin)){
		if(fscanf(fin, "%s", buffer) != 1){
			printf("Error: cannot read gene tree file 1 %s\n", inputfile);
			return ERROR;
		}
		
		for(char *p=buffer; *p; p++){
			*p=tolower(*p);
		} 

		if(!strcmp(buffer,"end;")){
			break;
		}

		if(!strcmp(buffer,"trees;") || !strcmp(buffer,"trees") || !strcmp(buffer,"tree")){
			n = 0;
		}else{
			n++;
		}
	}

	if (numGenes > n){
		printf("There are %d < %d gene trees in the input file %s\n",n, numGenes, inputfile);
		return ERROR;
	}

	fclose(fin);
	return NO_ERROR;
}

int ReadSPtableFromInputFile (char *inputfile){
	int index, i, j, n;
	char buffer[LSPNAME];	
	FILE *fin = fopen(inputfile, "r");

	if(!fin){
		printf("Error: cannot open file %s\n", inputfile);
		return ERROR;
	}

	while(!feof(fin)){
		if(fscanf(fin, "%s", buffer) != 1){
			printf("Error: cannot read input file %s\n", inputfile);
			return ERROR;
		}
		for(char *p=buffer; *p; p++){
			*p=tolower(*p);
		} 
		
		if(!strcmp(buffer,"matrix")){
			index = 0;
			for(i=0; i<sptree.ntaxa; i++){
				if(fscanf(fin,"%s %d", sptable.speciesname[i], &n) != 2){
					printf("Errors in sptable in file %s\n", inputfile);
					return ERROR;
				}
				for(j=0; j<n; j++){
					if(fscanf(fin,"%s", sptable.taxaname[index]) != 1){
						printf("Errors in sptable in file %s\n", inputfile);
						return ERROR;
					}	
					/*defines the species number that each taxon belongs to*/ 
					sptable.taxanodenumber[index++] = i;
				}
			}
			break;			
		}
	}
	fclose(fin);

	return NO_ERROR;
}

int InitialSpeciesTreePara(){
	int i;

	/*species names*/
	for(i=0; i<sptree.ntaxa; i++){
        strcpy(sptree.nodes[i].taxaname, sptable.speciesname[i]);
    }

	/*find outgroupnode*/
	for(i=0; i<sptree.ntaxa; i++){
		if(!strcmp(sptree.nodes[i].taxaname, sptree.outgroup)){
			break;
		}
	}
	if(i == sptree.ntaxa){
		printf("cannot find the outgroup species %s\n", sptree.outgroup);
		return ERROR;
	}
	sptree.outgroupnode = i;

	/*initialize other parameters*/
	for(i=0; i<2*sptree.ntaxa-1; i++){
		sptree.nodes[i].theta = 0.0;
		sptree.nodes[i].support = 0.0;
		if(i<sptree.ntaxa){
			sptree.nodes[i].brlens = 1.0;
		}
	}		

	return NO_ERROR;
}

int UpdateSPtable (){
	int i, j, k;

	/*update sptable taxanodenumber*/
	for(i=0; i<totaltaxa; i++){
		for(j=0; j<sptree.ntaxa; j++){
			if(sptable.taxanodenumber[i] == j){
				for(k=0;k<sptree.ntaxa;k++){
					if(strcmp(sptable.speciesname[j],sptree.nodes[k].taxaname)==0)  {
						sptable.taxanodenumber[i] = k;
						break;
					}
				}
				if(k == sptree.ntaxa){
					printf("%s is missing in the user tree\n\n",sptable.speciesname[j]);
					return ERROR;
				}
				break;
			}
		}
	}

	/*update sptable speciesname*/
	for(i=0; i<sptree.ntaxa; i++){
		strcpy(sptable.speciesname[i], sptree.nodes[i].taxaname);
	}

	/*update outgroupnode*/
	for(i=0; i<sptree.ntaxa; i++){
		if(!strcmp(sptree.nodes[i].taxaname, sptree.outgroup)){
			break;
		}
	}
	if(i == sptree.ntaxa){
		printf("cannot find the outgroup species %s\n", sptree.outgroup);
		return ERROR;
	}
	sptree.outgroupnode = i;

	return NO_ERROR;
}

int NJst (double *nodedistance, char **speciesname, uint32_t species_count, int outgroupSpecies) {
	int i, j, index;
	double *element;
	dist_matrix *dmat;
	
	dmat = dist_matrix_init(species_count);
    if (!dmat) {
        perror("Unable to create distance matrix");
        return NO_ERROR;
    }

	index = 0;
    for (i = 0; i < species_count; i++) {
		dist_matrix_set_species_name(dmat, i, speciesname[i]);
        dmat->cluster_sizes[i] = 1;
		if(i==outgroupSpecies){/*add a large distance for outgroupSpecies to root the tree*/
			for (j = 0; j < i; j++) {
				/*dist_matrix is the lower triangler matrix*/
				element = dist_matrix_element(dmat, i, j);
				*element = 2*species_count;
			}
		}else{
			for (j = 0; j < i; j++) {
				/*dist_matrix is the lower triangler matrix*/
				element = dist_matrix_element(dmat, i, j);
				*element = nodedistance[index++];
			}
		}
    }

	/*building the NJ tree using the function from neighbour_joining.c*/
	NeighborJoining(dmat);
	
	return NO_ERROR;	
}

int GenetreeNodeDistance (char *treefile, double *nodedistance){
	int numdistance = sptree.ntaxa * (sptree.ntaxa-1) / 2;
	int i, j, k, index, stop, sp1, sp2, location, *muldistcount;
	Tree genetree;
	FILE *fin;
	double *distance;
	
	/*allocate memory for distance; muldistcount for multiple alleles per species*/
	distance = (double*) malloc (numdistance * sizeof(double*));
	muldistcount = (int*) malloc (numdistance * sizeof(int*));
	for(i=0; i<numdistance; i++){
		distance[i] = 0.0;
		muldistcount[i] = 0;
	}

	fin = fopen(treefile,"r");
	if(!fin){
		printf("Error: cannot open file %s\n", treefile);
		return ERROR;
	}
	/*find node distances*/
	for(i=0; i<numGenes; i++){
		if(ReadaTree(fin, &genetree) == ERROR) {
			printf("Errors in gene tree %d.\n", i+1);
			return ERROR;
		}

		/*find gene tree namenumber*/
		for(j=0; j<genetree.ntaxa; j++){
			stop = 0;
			for(k=0; k<totaltaxa; k++){
				if(!strcmp(genetree.nodes[j].taxaname, sptable.taxaname[k])){
					stop = 1; 
					genetree.nodes[j].namenumber = sptable.taxanodenumber[k];
					break;
				}
			}
			if(stop == 0){
				printf("The taxon %s in gene tree %d is missing in the species-allele table\n", genetree.nodes[j].taxaname, i);
				return ERROR;
			}
		}

		for(j=0; j<genetree.ntaxa-1; j++){
			for(k=j+1; k<genetree.ntaxa; k++){
				sp1 = genetree.nodes[j].namenumber;
				sp2 = genetree.nodes[k].namenumber;

				/*calculate node distance for two distinct species*/
				if(sp1 != sp2){
					if(NodeDistanceLocation (sp1, sp2, &location) == ERROR){
						printf("Errors in NodeDistanceLocation\n");
						return ERROR;
					}
					distance[location] += NodeDistance(j, k, &genetree);
					muldistcount[location] += 1;
				}				
			}
		}
	}

	index=0;
	for(j=0; j<sptree.ntaxa-1; j++){
		for(k=j+1; k<sptree.ntaxa; k++){
			nodedistance[index] = distance[index]/muldistcount[index];
			if(nodedistance[index] >= sptree.ntaxa){
				printf("nodedistance is greater than sptree.ntaxa\n");
				return ERROR;
			}
			index++;
		}
	}

	free(distance);
	free(muldistcount);
	fclose(fin);
	return NO_ERROR;
}

int NodeDistanceLocation (int sp1, int sp2, int *location){ /*lower triangler matrix*/
	int node1, node2;

	if(sp1 < sp2){
		node1 = sp1;
		node2 = sp2;
	}else if (sp1 > sp2){
		node1 = sp2;
		node2 = sp1;
	}else{
		printf("Errors in NodeDistanceLocation: two identical species!\n");
		return ERROR;
	}	
	*location = node2*(node2-1)/2+node1;

	return NO_ERROR;
}

int GenetreePartitions (char *treefile, int ntrees){
	int numpartitions = sptree.ntaxa-1;
	int i, j, k, x, stop;
	Tree genetree;
	FILE *fin, *fout;
	char outfile[100];
	int **genetreepartitions, **genetreepartitions1;

	
	fin = fopen(treefile,"r");
	if(!fin){
		printf("Error: cannot open file %s\n", treefile);
		return ERROR;
	}

	sprintf(outfile, "%s.partition", treefile);
	fout = fopen(outfile,"w");
	if(!fout){
		printf("Error: cannot open file %s\n", outfile);
		return ERROR;
	}

	/*allocate memory for genetreepartitions*/
	genetreepartitions = (int**)calloc(numpartitions, sizeof(int*));
    genetreepartitions[0] = (int*)calloc(sptree.ntaxa * numpartitions, sizeof(int));
   	for(i = 0; i < numpartitions; i++){
		genetreepartitions[i] = genetreepartitions[0] + i * sptree.ntaxa;
	}
  	if(!genetreepartitions){
		printf("allocating problems for genetreepartitions!\n");
	   	return ERROR;
	} 

	/*allocate memory for genetreepartitions*/
	genetreepartitions1 = (int**)calloc(numpartitions, sizeof(int*));
    genetreepartitions1[0] = (int*)calloc(sptree.ntaxa * numpartitions, sizeof(int));
   	for(i = 0; i < numpartitions; i++){
		genetreepartitions1[i] = genetreepartitions1[0] + i * sptree.ntaxa;
	}
  	if(!genetreepartitions1){
		printf("allocating problems for genetreepartitions1!\n");
	   	return ERROR;
	} 

	/*find gene tree partitions*/
	for(i=0; i<ntrees; i++){
		if(ReadaTree(fin, &genetree) == ERROR) {
			printf("Errors in the gene tree %d.\n", i+1);
			return ERROR;
		}

		/*find gene tree namenumber*/
		for(j=0; j<genetree.ntaxa; j++){
			stop = 0;
			for(k=0; k<totaltaxa; k++){
				if(!strcmp(genetree.nodes[j].taxaname, sptable.taxaname[k])){
					stop = 1; 
					genetree.nodes[j].namenumber = sptable.taxanodenumber[k];
					break;
				}
			}
			if(stop == 0){
				printf("The taxon %s in gene tree %d is missing in the species-allele table in the control file\n", genetree.nodes[j].taxaname, i);
				return ERROR;
			}
		}

		/*refresh partition matrix*/
		for(j=0; j<numpartitions; j++){
			for(k=0; k<sptree.ntaxa; k++){
				genetreepartitions[j][k] = 0;
				genetreepartitions1[j][k] = 0;
			}
		}

		if(TreePartitions (genetreepartitions, &genetree) == ERROR){
			printf("Errors in the TreePartitions function\n");
			return ERROR;
		}

		for(j=0; j<numpartitions; j++){
			for(k=0; k<sptree.ntaxa; k++){
				if(genetreepartitions[j][k] == 1){
					x = genetree.nodes[k].namenumber;
					genetreepartitions1[j][x] = 1;
				}
			}
		}

		fprintf(fout,"tree%d\t",i+1);
		for(j=0; j<numpartitions; j++){
			for(k=0; k<sptree.ntaxa; k++){
				fprintf(fout,"%d", genetreepartitions[j][k]);
			}
			fprintf(fout," ");			
		}
		fprintf(fout,"\n");
	}

	free(genetreepartitions[0]);
	free(genetreepartitions);
	free(genetreepartitions1[0]);
	free(genetreepartitions1);
	fclose(fin);
	fclose(fout);
	
	return NO_ERROR;
}

int SptreePartitionSupport (char *treefile){
	int numpartitions = sptree.ntaxa-1;
	int i, j, k, l, x, stop;
	Tree genetree;
	FILE *fin;
	int **genetreepartitions, **genetreepartitions1, **partitions;
	
	/*allocate memory for partitions*/
	partitions = (int**)calloc(numpartitions, sizeof(int*));
    partitions[0] = (int*)calloc(sptree.ntaxa * numpartitions, sizeof(int));
   	for(i = 0; i < numpartitions; i++){
		partitions[i] = partitions[0] + i * sptree.ntaxa;
	}
  	if(!partitions){
		printf("allocating problems for partitions!\n");
	   	return ERROR;
	} 

	/*allocate memory for genetreepartitions*/
	genetreepartitions = (int**)calloc(numpartitions, sizeof(int*));
    genetreepartitions[0] = (int*)calloc(sptree.ntaxa * numpartitions, sizeof(int));
   	for(i = 0; i < numpartitions; i++){
		genetreepartitions[i] = genetreepartitions[0] + i * sptree.ntaxa;
	}
  	if(!genetreepartitions){
		printf("allocating problems for genetreepartitions!\n");
	   	return ERROR;
	} 

	/*allocate memory for genetreepartitions*/
	genetreepartitions1 = (int**)calloc(numpartitions, sizeof(int*));
    genetreepartitions1[0] = (int*)calloc(sptree.ntaxa * numpartitions, sizeof(int));
   	for(i = 0; i < numpartitions; i++){
		genetreepartitions1[i] = genetreepartitions1[0] + i * sptree.ntaxa;
	}
  	if(!genetreepartitions1){
		printf("allocating problems for genetreepartitions1!\n");
	   	return ERROR;
	} 

	fin = fopen(treefile,"r");
	if(!fin){
		printf("Error: cannot open file %s\n", treefile);
		return ERROR;
	}

	/*find species tree partitions*/
	TreePartitions(partitions, &sptree);

	/*refresh species tree support*/
	for(j=sptree.ntaxa; j<2*sptree.ntaxa-1; j++){
		sptree.nodes[j].support = 0.0;
	}

	/*find gene tree partitions*/
	for(i=0; i<numGenes; i++){
		if(ReadaTree(fin, &genetree) == ERROR) {
			printf("Errors in the gene tree %d.\n", i+1);
			return ERROR;
		}

		/*find gene tree namenumber*/
		for(j=0; j<genetree.ntaxa; j++){
			stop = 0;
			for(k=0; k<totaltaxa; k++){
				if(!strcmp(genetree.nodes[j].taxaname, sptable.taxaname[k])){
					stop = 1; 
					genetree.nodes[j].namenumber = sptable.taxanodenumber[k];
					break;
				}
			}
			if(stop == 0){
				printf("The taxon %s in gene tree %d is missing in the species-allele table in the control file\n", genetree.nodes[j].taxaname, i);
				return ERROR;
			}
		}

		/*refresh partition matrix*/
		for(j=0; j<numpartitions; j++){
			for(k=0; k<sptree.ntaxa; k++){
				genetreepartitions[j][k] = 0;
				genetreepartitions1[j][k] = 0;
			}
		}

		if(TreePartitions (genetreepartitions, &genetree) == ERROR){
			printf("Errors in the TreePartitions function\n");
			return ERROR;
		}

		/*change to species namenumber*/
		for(j=0; j<numpartitions; j++){
			for(k=0; k<sptree.ntaxa; k++){
				if(genetreepartitions[j][k] == 1){
					x = genetree.nodes[k].namenumber;
					genetreepartitions1[j][x] = 1;
				}
			}
		}

		/*compare genetreepartitions1 with species tree partitions*/
		for(j=0; j<numpartitions; j++){
			for(k=0; k<numpartitions; k++){
				x = 0;
				if(partitions[j][0]==genetreepartitions1[k][0]){
					for(l=1; l<sptree.ntaxa; l++){
						if(partitions[j][l] != genetreepartitions1[k][l]){
							x = 1;
							break;
						}
					}
				}else{
					for(l=1; l<sptree.ntaxa; l++){
						if(partitions[j][l] == genetreepartitions1[k][l]){
							x = 1;
							break;
						}
					}
				}
				if(x==0){					
					sptree.nodes[j+sptree.ntaxa].support += 1.0/numGenes;
					break;
				}
			}				
		}	
	}

	free(genetreepartitions[0]);
	free(genetreepartitions);
	free(genetreepartitions1[0]);
	free(genetreepartitions1);
	fclose(fin);
	
	return NO_ERROR;
}


/*postorder: right, left, root*/
int	TreeTraversalPostorder (int node, Tree *tree){
	int son, i;
	
	if(node < tree->ntaxa){
		tree->postorder[postorderindex++] = node;
	}else{
		/*travel in the reverse order*/
		for(i=tree->nodes[node].nson-1; i>=0; i--){
			son = tree->nodes[node].sons[i];
			TreeTraversalPostorder(son, tree);
		}
		tree->postorder[postorderindex++] = node;
	}
	return NO_ERROR;
}

int TreePartitions (int **partitions, Tree *tree){
	int i, index;
	
	index = 0;
	for(i=tree->ntaxa; i<tree->numnodes; i++){
		FindOffsprings(partitions[index++], tree, i);
	}

	return NO_ERROR;
}

/*preorder: root, left, right*/
int	TreeTraversalPreorder (int node, Tree *tree){
	int son, i;
	if(node < tree->ntaxa){
		tree->preorder[postorderindex++] = node;
	}else{
		tree->preorder[postorderindex++] = node;
		for(i=0; i<tree->nodes[node].nson; i++){
			son = tree->nodes[node].sons[i];
			TreeTraversalPreorder(son, tree);
		}	
	}
	return NO_ERROR;
}

int TreeBranchCollapse (char *treefile){
   	int i, j, index, x, nleft=0, nright=0;
	int leftposition[MAX_STRING_LENGTH], rightposition[MAX_STRING_LENGTH], position[MAX_STRING_LENGTH], treestrlength;
   	char treestring[MAX_STRING_LENGTH], str[100], ch, cha=' ';
	double brlens;
	char outfile[100];
	FILE *fin, *fout;

	fin = fopen(treefile, "r");
	if(!fin){
		printf("Error: cannot open file %s\n", treefile);
		return ERROR;
	}

	sprintf(outfile, "%s_collapse.tre", treefile);
	fout = fopen(outfile,"w");
	if(!fout){
		printf("Error: cannot open file %s\n", outfile);
		return ERROR;
	}

	while(cha != EOF){
		/*move to the first character of a treestring*/
		while(cha != '(' && cha != EOF){
			cha = fgetc(fin);
		}

		if(cha == EOF){
			break;
		}

		/*read treestring*/	
		index = 0;
		while(cha != ';'){		
			treestring[index++] = cha;
			cha = fgetc(fin);

			/*skip white space*/
			while(cha == ' '){
				cha = fgetc(fin);
			}	
		}

		/*add ending to treestring*/
		treestring[index++] = ';';
		treestrlength = index;
		treestring[treestrlength] = '\0';

		nleft = nright = 0;
		for(i=0; i<treestrlength; i++){
			if(treestring[i] == '('){
				leftposition[nleft++] = i;
			}else if (treestring[i] == ')'){
				rightposition[nright++] = i;
			}else{
				continue;
			}
		}

		/*check if the number of left == the number of right parentheses*/
		if(nleft != nright){
			printf("Errors in TreeBranchCollapse");
			return ERROR;
		}


		/*for each right, we find the corresponding left parenthesis stored in position*/
		for(i=0; i<nright; i++){
			index = 1;
			while(leftposition[nleft-index] > rightposition[i]){
				index++;
			}
			position[i] = leftposition[nleft-index];
			leftposition[nleft-index] = rightposition[nright-1]+1;
		}	

		for(i=0; i<nright; i++){
			x = rightposition[i];
			if(treestring[x+1] == ':'){
				if(!(isdigit(treestring[x+2]) || treestring[x+2] == '.')){
					printf("Errors in TreeBranchCollapse\n");
					return ERROR;
				}

				ch = treestring[++x];

				index = 0;
				while(ch != ',' && ch != ')' && ch != ';'){
					ch = treestring[++x];
					str[index++] = ch;	
				}
				x--;
				str[--index] = '\0';

				//printf("str %s %d %d %c\n", str, x, index, treestring[x]);
				
				sscanf(str, "%lf", &brlens);
				if(brlens != NA && brlens < collapseBrlens){
					for(j=x-index; j<x+1; j++){
						treestring[j] = '*';
					}
					treestring[rightposition[i]] = treestring[position[i]] = '*';
				}
			}	
		}
		for(i=0; i<treestrlength; i++){
			if(treestring[i] != '*'){
				fprintf(fout, "%c", treestring[i]);
			}
		}
		fprintf(fout, "\n");
	}
	
	fclose(fin);
	fclose(fout);
	return NO_ERROR;
}

int Loglikelihood (int **triple, Tree *tree, double *loglike){
	int i, j, k, w, son0, son1, father, cnode, ntriples=0;
	int offsprings0[NTAXA], offsprings1[NTAXA], offsprings2[NTAXA], taxa0[NTAXA], taxa1[NTAXA], taxa2[NTAXA], n0, n1, n2;
	int location[2];
	double p, brlens, logp;
	int *array;
		
	*loglike = 0.0;
	for(i=tree->ntaxa; i<2*tree->ntaxa-1; i++){
		if(i != tree->root){
			brlens = tree->nodes[i].brlens;
			n0 = n1 = 0;
			for(j=0; j<tree->ntaxa; j++){
				offsprings0[j] = 0;
				offsprings1[j] = 0;
			}

			son0 = tree->nodes[i].sons[0];
			son1 = tree->nodes[i].sons[1];
			FindOffsprings (offsprings0, tree, son0);
			FindOffsprings (offsprings1, tree, son1);

			for(j=0; j<tree->ntaxa; j++){
				if (offsprings0[j] == 1)
					taxa0[n0++] = j;
				if (offsprings1[j] == 1)
					taxa1[n1++] = j;
			}

			cnode = i;
			father = tree->nodes[i].father;
			while (father != -1){
				for(j=0; j<tree->ntaxa; j++){
					offsprings2[j] = 0;
				}
											
				n2 = 0;
				p = 1-exp(-brlens)*2/3;

				if(tree->nodes[father].sons[0] == cnode){
					FindOffsprings (offsprings2, tree, tree->nodes[father].sons[1]);
				}else{
					FindOffsprings (offsprings2, tree, tree->nodes[father].sons[0]);
				}

				for(j=0; j<tree->ntaxa; j++){
					if (offsprings2[j] == 1){
						taxa2[n2++] = j;
					}		
				}

				for(j=0; j<n0; j++){
					for(k=0; k<n1; k++){
						for(w=0; w<n2; w++){
							FindTriple (taxa0[j], taxa1[k], taxa2[w], location);
							if (DEBUG){
								array[location[0]] = 1;
								/*printf("location %d %d %d %ld %ld %f %d %d\n",taxa0[j],taxa1[k],taxa2[w],location[0],location[1],brlens,triplematrix[location[0]][3],triplematrix[location[0]][location[1]]);*/
							}
							if(LogBinomialP (triplematrix[location[0]][3],triplematrix[location[0]][location[1]],p,&logp) == ERROR){
								printf("Errors in LogBinomialP\n");
								return ERROR;
							}
							/*printf("logp %f %d %d %d\n",logp, taxa0[j],taxa1[k],taxa2[w]);*/
							*loglike += logp;
						}
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

	if(ntriples != (tree->ntaxa * (tree->ntaxa-1) * (tree->ntaxa-2) / 6)){
		return ERROR;
	}

	return NO_ERROR;
}

int LogBinomialP (int n, int x, double p, double *logp){
	if(p == 1) p = 0.99999999;
	*logp = x*log(p) + (n-x)*log((1-p)/2);
	return NO_ERROR;
}

int Algorithm (int **triple, FILE *besttreefile, FILE *outputfile, int runindex){
	double oldloglike, diff;
	int i, round=0, sample = 1000, totalround = MAXROUND, noincreasing=0, accept = 0;
	Tree oldsptree;
	time_t starttime, endtime;
	clock_t start, end;
    double cpu_time_used;
	
    /*starting tree*/
	if(usertree == 0){
		//RandomTree (&sptree);
	} 

	CopyTree (&sptree, &oldsptree);
	for(i=0; i<sptree.ntaxa; i++){
		strcpy(oldsptree.nodes[i].taxaname, sptree.nodes[i].taxaname);
	} 
	
	if(Loglikelihood (triple, &sptree, &curLn) == ERROR){
		printf("Errors in Loglikelihood\n");
		return ERROR;
	}

    /*calculating the likelihood score for a fixed tree, which is given as the user tree*/
    if(usertree == 2) {
		printf("loglikelihood of the species tree is %lf\n", curLn);
		fprintf(outputfile, "Begin loglikelihood;\n  loglike = %lf;\n", curLn);
		return NO_ERROR;
	}

    /*record cpu time*/
	start = clock();
	starttime = time(0);

	while (round < totalround){
		/*copy likelihood and tree*/
		CopyTree (&sptree, &oldsptree);
		oldloglike = curLn;
		
		if(MoveNode (&sptree) == ERROR){
			printf("Errors in MoveNode\n");
			return ERROR;
		}
		if(Loglikelihood (triple, &sptree, &curLn) == ERROR){
			printf("Errors in Loglikelihood\n");
			return ERROR;
		}
		diff = curLn - oldloglike;

		if(diff < 0){
			noincreasing++;
			curLn = oldloglike;
			CopyTree (&oldsptree, &sptree);
		}else{
			if(diff > 0)
				noincreasing = 0;			
			accept++;
		}
		
		round ++;

		/*screen printout*/
		if((round % sample) == 0 || round == 1){
			printf("\tround %d \t\tloglike %f \t....completed....\n", round, curLn);
		}

		if(runindex == 0){
			/*print out all trees generated in the algorithm. if nrun > 1, we only print out trees for the first run*/
			if((round % sample) == 0 || round == 1){
				if (PrintState (round, outputfile, 0) == ERROR){
					Print ("%S   Errors in PrintState.\n", spacer);
					return ERROR;
				}
			}

			/*print out the final mpest tree*/
			if(runindex == 0 && round == 1){
				if (PrintState (round, besttreefile, 1) == ERROR){
					Print ("%S   Errors in PrintState.\n", spacer);
					return ERROR;
				}
			}
		}
			
		if (noincreasing > NUM_NOCHANGE || noincreasing == NUM_NOCHANGE){
			/*screen output*/
			printf("\tround %d \t\tloglike %f \t....completed....\n", round, curLn);

			/*if nrun > 1, we only print out trees for the first run*/
			if(runindex == 0){
				if (PrintState (round, outputfile, 0) == ERROR){
					Print ("%S   Errors in PrintState.\n", spacer);
					return ERROR;
				}
			}

			/*print out the final mpest tree for all runs*/
			if (PrintState (round, besttreefile, 1) == ERROR){
				Print ("%S   Errors in PrintState.\n", spacer);
				return ERROR;
			}
			break;
		} 	
	}
	
	/*record cpu time*/
	end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	endtime = time(0);

	if (difftime (endtime, starttime) > 2.0){
		printf("\t\t\tAnalysis is completed in %.0f seconds (cpu time is %.0f seconds)\n", difftime(endtime, starttime), cpu_time_used);
	}else if (difftime (endtime, starttime) >= 1.0){
		printf("\t\t\tAnalysis completed in 1 second\n");
	}else{
		printf("\t\t\tAnalysis completed in less than 1 second\n");
	}

	return NO_ERROR;
}

int PrintState (int round, FILE *outfile, int addend){
	int i;
	char buffer[30];
  	struct timeval tv;
  	time_t curtime; 

	/*print to file*/
	if (addend == 0){
		if(round == 1){
			gettimeofday(&tv, NULL);
			curtime = tv.tv_sec;
			strftime(buffer,30,"%I:%M%p on %m-%d-%Y",localtime(&curtime));

			fprintf(outfile, "#Nexus\n[This analysis was conducted at local time %s with seed = %ld.\nThe branch length 1.0 of the external branches is arbitrary because mp-est does not\nestimate the lengths of the external branches. If the internal branch length is 7.0,\nit indicates that all gene trees support the species tree triple. In the output.tre file,\nthe tree block contains the output trees (and their likelihood scores) generated from the algorithm,\nwhile the besttree.tre file contains the final mpest tree. If multiple runs are specified,\nthe besttree.tre file contains multiple mpest trees generated from multiple runs.\nChoose the tree with the maximum likelihood score as the estimate of the species tree.]\n\n", buffer, seed);
			fprintf(outfile, "Begin trees;\n  translate\n");
			for (i=0; i<sptree.ntaxa-1; i++){
				fprintf(outfile,"    %d %s,\n", i+1, sptree.nodes[i].taxaname);
			}
			fprintf(outfile,"    %d %s;\n", sptree.ntaxa, sptree.nodes[sptree.ntaxa-1].taxaname);
			if (PrintTree(&sptree, sptree.root, 0, 1, 0, 0, 1) == ERROR){
                printf("Errors in printtree!\n");
                return ERROR;
            }
            fprintf(outfile, "  tree round%d [%2.6f] = %s", round, curLn, printString);
            free (printString);
        }else{
            if (PrintTree(&sptree, sptree.root, 0, 1, 0, 0, 1) == ERROR){
                printf("Errors in printtree!\n");
                return ERROR;
            }
            fprintf(outfile, "  tree round%d [%2.6f] = %s", round, curLn, printString);
            free (printString);
        }
	}else{
		if(round == 1){
			gettimeofday(&tv, NULL);
			curtime = tv.tv_sec;
			strftime(buffer,30,"%I:%M%p on %m-%d-%Y",localtime(&curtime));

			fprintf(outfile, "#Nexus\n[This analysis was conducted at local time %s with seed = %ld.\nThe branch length 1.0 of the external branches is arbitrary because mp-est does not\nestimate the lengths of the external branches. If the internal branch length is 7.0,\nit indicates that all gene trees support the species tree triple. In the output.tre file,\nthe tree block contains the output trees (and their likelihood scores) generated from the algorithm,\nwhile the besttree.tre file contains the final mpest tree. If multiple runs are specified,\nthe besttree.tre file contains multiple mpest trees generated from multiple runs.\nChoose the tree with the maximum likelihood score as the estimate of the species tree.]\n\n", buffer, seed);
			fprintf(outfile, "Begin trees;\n  translate\n");
			for (i=0; i<sptree.ntaxa-1; i++) {
				fprintf(outfile,"    %d %s,\n", i+1, sptree.nodes[i].taxaname);
			}
			fprintf(outfile,"    %d %s;\n", sptree.ntaxa, sptree.nodes[sptree.ntaxa-1].taxaname);
        }else{
			/*calculate gene tree frequencies for 1 allele per species*/
			if(totaltaxa == sptree.ntaxa){
				SptreePartitionSupport (genetreefile);
			}
			if (PrintTree(&sptree, sptree.root, 0, 1, 0, 1, 1) == ERROR){
				printf("Errors in printtree!\n");
				return ERROR;
			}
			fprintf(outfile, "  tree mpest [%2.6f] = %s", curLn, printString);
			free (printString);
		}
	}
	return NO_ERROR;
}

/*this function is for debugging purpose*/
int printSptree (void){
	int i;

    if (PrintTree(&sptree, sptree.root, 1, 1, 0, 0, 1) == ERROR){
		printf("Errors in printtree!\n");
		return ERROR;
	}
	printf("%s", printString);
	free (printString);

	for(i=0; i<2*sptree.ntaxa-1; i++){
		printf("node %d %s %d %d %d\n", i, sptree.nodes[i].taxaname, sptree.nodes[i].father,sptree.nodes[i].sons[0],sptree.nodes[i].sons[1]);
	}
		
    return NO_ERROR;
}

int MoveNode (Tree *tree){
	int inode, jnode, father, grandfather, nodes[5], nnodes=0, son;
	double rand = rndu();
	double p[2]={0.3,0.8};

	if(usertree == 3) {p[0] = 0.0; p[1]=0.0;}
    if(usertree == 4) {p[0] = 0.0; p[1] = 0.7;}
    
	if(rand < p[0]){
		do{
			do{
				inode = (int)(rndu() * (2*tree->ntaxa-1));
			}while (inode == tree->root || inode == tree->nodes[tree->root].sons[0] || inode == tree->nodes[tree->root].sons[1]);

			father = tree->nodes[inode].father;
			grandfather = tree->nodes[father].father;

			if(tree->nodes[father].sons[0] == inode){
				son = tree->nodes[father].sons[1];
			}else{
				son = tree->nodes[father].sons[0];
			}
				

			if(son >= tree->ntaxa){
				nodes[nnodes++] = tree->nodes[son].sons[0];
				nodes[nnodes++] = tree->nodes[son].sons[1];
			}		
				
			if(tree->nodes[grandfather].sons[0] == father){
				nodes[nnodes++] = son = tree->nodes[grandfather].sons[1];
			}else{
				nodes[nnodes++] = son = tree->nodes[grandfather].sons[0];
			}
		
			if(son >= tree->ntaxa){
				nodes[nnodes++] = tree->nodes[son].sons[0];
				nodes[nnodes++] = tree->nodes[son].sons[1];
			}
	    }while(nnodes == 0);

		jnode = nodes[(int)(rndu() * nnodes)];
	
		SwapNodes (tree, inode, jnode);
	}else if (rand < p[1]){
		do{
			inode = rndu() * tree->ntaxa;
            if(usertree == 4) inode = updatenodes[(int)(rndu()*numupdatenodes)];
		}while (tree->nodes[inode].father == tree->root);

        /*delete the branch of inode*/
	    if(DeleteNode(tree, inode) == ERROR){
		 printf("Errors in DeleteNode\n");
		 return ERROR;
	    }

		do{
			jnode = rndu() * (2*tree->ntaxa-1);
		}while (jnode == tree->root || jnode == tree->nodes[inode].father || jnode == inode);

        /*attach inode to the branch of jnode*/
        if(AddNode(tree, inode, tree->nodes[inode].father, jnode) == ERROR){
		  printf("Errors in AddNode\n");
		  return ERROR;
	    }
	}else{
		MoveBrlens (tree);
	}

	return NO_ERROR;
}

void MoveBrlens (Tree *tree){
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

void CopyTree(Tree *from, Tree *to){
   	int   i ;
 
  	/*copy the species tree nodes*/
  	to->root = from->root;
	to->ntaxa = from->ntaxa;
  	for(i=0; i<(2*from->ntaxa-1); i++){
    		to->nodes[i].nson = from->nodes[i].nson;
    		to->nodes[i].sons[0] = from->nodes[i].sons[0];
    		to->nodes[i].sons[1] = from->nodes[i].sons[1];
    		to->nodes[i].father = from->nodes[i].father;
    		to->nodes[i].brlens = from->nodes[i].brlens;
  	}
}

int FindOutgroup (Tree *tree){
	int outgroup;

	if(tree->nodes[tree->root].sons[0] < tree->ntaxa)
		outgroup = tree->nodes[tree->root].sons[0];
	else
		outgroup = tree->nodes[tree->root].sons[1];
	return (outgroup);
}

int DeleteNode (Tree *tree, int inode){
	int father, son, grandfather;
	
	if(inode == tree->root) return ERROR;

	father = tree->nodes[inode].father;

	if(tree->nodes[father].sons[0] == inode){
		son = tree->nodes[father].sons[1];
	}else{
		son = tree->nodes[father].sons[0];
	} 

	if(father == tree->root){
		tree->root = son;
	}else{
		grandfather = tree->nodes[father].father;
		if(tree->nodes[grandfather].sons[0] == father)
			tree->nodes[grandfather].sons[0] = son;
		else
			tree->nodes[grandfather].sons[1] = son;
		tree->nodes[son].father = grandfather;
	}

	return NO_ERROR;
}

int AddNode (Tree *tree, int fromnode, int fromfather, int tonode){
	int tofather;

    /*fromnode should not be the outgroup species*/
	if(fromfather == -1){
		printf("Errors for fromfather\n");
		return ERROR;
	}

    /*update the ancestral informatino of fromfather*/
	tofather = tree->nodes[tonode].father;
	if(tree->nodes[tofather].sons[0] == tonode){
		tree->nodes[tofather].sons[0] = fromfather;
	}else{
		tree->nodes[tofather].sons[1] = fromfather;
	}	
	tree->nodes[fromfather].father = tofather;

    /*update the ancestral information of tonode*/
	tree->nodes[tonode].father = fromfather;
	if(tree->nodes[fromfather].sons[0] == fromnode){
		tree->nodes[fromfather].sons[1] = tonode;
	}else{
		tree->nodes[fromfather].sons[0] = tonode;
	} 

	return NO_ERROR;
}

int SwapNodes (Tree *tree, int inode, int jnode){
	int ifather, jfather;

	ifather = tree->nodes[inode].father;
	jfather = tree->nodes[jnode].father;
	tree->nodes[inode].father = jfather;
	tree->nodes[jnode].father = ifather;

	if(tree->nodes[ifather].sons[0] == inode){
		tree->nodes[ifather].sons[0] = jnode;
	}else{
		tree->nodes[ifather].sons[1] = jnode;
	}
		
	if(tree->nodes[jfather].sons[0] == jnode){
		tree->nodes[jfather].sons[0] = inode;
	}else{
		tree->nodes[jfather].sons[1] = inode;
	}
	
	return NO_ERROR;
}

void RandomVector (int *array, int number){
	int i, n, m;

	for(i=0; i<number; i++){
		array[i] = i;
	}

	for(i=0; i<number; i++){
		n = rndu() * (number - i);
		m = array[number-i-1];
		array[number-i-1] = array[n];
		array[n] = m;		
	}
}

void RandomTree (Tree *tree){
	int i, outgroup, array[NTAXA];
	double internalbrlens = 0.1, externalbrlens = 1.0;

	outgroup = 0;
	RandomVector (array, tree->ntaxa);
	tree->nodes[array[0]].nson = 0;
	tree->nodes[array[0]].sons[0] = -1;
	tree->nodes[array[0]].sons[1] = -1;
	tree->nodes[array[0]].father = tree->ntaxa;
	tree->nodes[tree->ntaxa].sons[0] = array[0];
	tree->nodes[array[0]].brlens = externalbrlens;
	
	tree->nodes[array[1]].nson = 0;
	tree->nodes[array[1]].sons[0] = -1;
	tree->nodes[array[1]].sons[1] = -1;
	tree->nodes[array[1]].father = tree->ntaxa;
	tree->nodes[tree->ntaxa].sons[1] = array[1];
	tree->nodes[array[1]].brlens = externalbrlens;
	
	tree->nodes[tree->ntaxa].father = tree->ntaxa + 1;
	tree->nodes[tree->ntaxa].brlens = internalbrlens;
	tree->nodes[tree->ntaxa].nson = 2;
	
	tree->nodes[2*tree->ntaxa-2].nson = 2;
	tree->nodes[2*tree->ntaxa-2].sons[1] = 2*tree->ntaxa-3;
	tree->nodes[2*tree->ntaxa-2].father = -1;
	tree->nodes[2*tree->ntaxa-2].brlens = 0;
	
	for(i=2; i<tree->ntaxa; i++){
		tree->nodes[array[i]].nson = 0;
		tree->nodes[array[i]].sons[0] = -1;
		tree->nodes[array[i]].sons[1] = -1;
		tree->nodes[array[i]].father = tree->ntaxa+i-1;
		tree->nodes[tree->ntaxa+i-1].sons[0] = array[i];
		tree->nodes[array[i]].brlens = externalbrlens;
	}

	for(i=tree->ntaxa+1; i<2*tree->ntaxa-2; i++){
		tree->nodes[i].nson = 2;
		tree->nodes[i].sons[1] = i-1;
		tree->nodes[i].father = i+1;
		tree->nodes[i].brlens = internalbrlens;
	}
	tree->root = 2*tree->ntaxa-2;
	tree->numnodes = 2*tree->ntaxa-1;
	tree->isrooted = 1;
	if(outgroup != array[tree->ntaxa - 1]){
		SwapNodes (tree, outgroup, array[tree->ntaxa-1]);
	}		
}

int TriplesFreq (char *treefile, int ngene){
	int i, j, k, stop;
	int ntriples = sptree.ntaxa*(sptree.ntaxa-1)*(sptree.ntaxa-2)/6;
	Tree genetree;
	FILE *fin;

	fin = fopen(treefile, "r");
	if(!fin){
		printf("Error: cannot open file %s\n", treefile);
		return ERROR;
	}
	
	/*find gene tree triples*/
	for(i=0; i<ngene; i++){
		if(ReadaTree(fin, &genetree) == ERROR) {
			printf("Errors in the gene tree %d; It must be a rooted tree.\n",i+1);
			return ERROR;
		}

		/*find gene tree namenumber*/
		for(j=0; j<genetree.ntaxa; j++){
			stop = 0;
			for(k=0; k<totaltaxa; k++){
				if(!strcmp(genetree.nodes[j].taxaname, sptable.taxaname[k])){
					stop = 1; 
					genetree.nodes[j].namenumber = sptable.taxanodenumber[k];
					break;
				}
			}
			if(stop == 0){
				printf("The taxon %s in gene tree %d is missing in the species-allele table in the control file\n", genetree.nodes[j].taxaname, i);
				return ERROR;
			}
		}

		if(TriplesFreqInatree (triplematrix, &genetree) == ERROR){
			printf("Errors in the TriplesFreqInatree function\n");
			return ERROR;
		}
	}

	for(i=0; i<ntriples; i++){
		for(j=0; j<3; j++){
			triplematrix[i][3] += triplematrix[i][j];
		}		
	}

	/*close file*/
	fclose(fin);

	return NO_ERROR;
}

void PrintHeader (void){
	Print ("\n\n\n\n%s            Maximum Pseudo-likelihood Estimation of Species Trees (MPEST_v2.1)  \n\n", spacer);
	Print ("%s                                   by\n\n", spacer);
	Print ("%s                                Liang Liu\n\n", spacer);
	Print ("%s                        Department of Statistics\n", spacer);
	Print ("%s                          University of Georgia\n", spacer);
	Print ("%s                               lliu@uga.edu\n\n", spacer);
	Print ("%s               Distributed under the GNU General Public License\n\n", spacer);	
}

int TriplesFreqInatree (int **triple, Tree *tree){
	int i, j, k, jj, kk, ww, son0, son1;
	int offsprings0[NTAXA], offsprings1[NTAXA], offsprings2[NTAXA], n0, n1, n2;
	int location[2];
 	
	for(i=tree->ntaxa; i<tree->numnodes; i++){
		if(i != tree->root){			
			for(j=0; j<tree->ntaxa; j++){
				offsprings2[j] = 0;
			}

			/*find the species not under node i*/
			FindOffsprings (offsprings2, tree, i);
			n2 = 0;
			for(j=0; j<tree->ntaxa; j++){
				if(offsprings2[j] == 0){
					offsprings2[n2++] = tree->nodes[j].namenumber;
				}
			}

			/*find the descendant species for son0 and son1*/
			for(j=0; j<tree->nodes[i].nson-1; j++){
				for(k=j+1; k<tree->nodes[i].nson; k++){
					son0 = tree->nodes[i].sons[j];
					son1 = tree->nodes[i].sons[k];

					/*reset offsprings0 and offsprings1*/
					for(jj=0; jj<tree->ntaxa; jj++){
						offsprings0[jj] = 0;
						offsprings1[jj] = 0;
					}

					FindOffsprings (offsprings0, tree, son0);
					FindOffsprings (offsprings1, tree, son1);

					/*find the species under each son group*/
					n0 = n1 = 0;
					for(jj=0; jj<tree->ntaxa; jj++){
						if(offsprings0[jj] == 0 && offsprings1[jj] == 0){
							continue;
							//offsprings2[n2++] = tree->nodes[j].namenumber;
						}else if (offsprings0[jj] == 1 && offsprings1[jj] == 0){
							offsprings0[n0++] = tree->nodes[jj].namenumber;
						}else if (offsprings0[jj] == 0 && offsprings1[jj] == 1){
							offsprings1[n1++] = tree->nodes[jj].namenumber;
						}else{
							printf("Error in TriplesFreqInatree\n");
							return ERROR;
						}
					}

					for(jj=0; jj<n0; jj++){
						for(kk=0; kk<n1; kk++){
							for(ww=0; ww<n2; ww++){
								/*triplets involving lineages from the same species are ignored. This occurs only when multiple alleles are sampled per species*/
								if(offsprings0[jj] != offsprings1[kk] && offsprings0[jj] != offsprings2[ww] && offsprings1[kk] != offsprings2[ww]){													
									FindTriple (offsprings0[jj], offsprings1[kk], offsprings2[ww], location);
									triple[location[0]][location[1]]++;
								}	
							}
						}
					}	
				}
			}					
		}
	}

	return NO_ERROR;
}

int FindNodeAncestors (int node, int *ancestors, Tree *tree){
	int numanc=1, k=node;
	ancestors[0] = node;
	while(k != tree->root){
		k = tree->nodes[k].father;
		ancestors[numanc++] = k;
	}
	return numanc;
}

int QuartetsFreqInatree (int **quartet, Tree *tree){
	int i, j, k, l, x1, x2, x3, x4;
	int dis1, dis2, dis3, dis4, dis5, dis6;
	int location[2];
 	
	for(i=0; i<tree->ntaxa-3; i++){
		for(j=i+1; j<tree->ntaxa-2; j++){
			for(k=j+1; k<tree->ntaxa-1; k++){
				for(l=k+1; l<tree->ntaxa; l++){
					x1 = tree->nodes[i].namenumber;
					x2 = tree->nodes[j].namenumber;
					x3 = tree->nodes[k].namenumber;
					x4 = tree->nodes[l].namenumber;
					
					if((x1-x2)*(x1-x3)*(x1-x4)*(x2-x3)*(x2-x4)*(x3-x4) != 0){
						dis1 = NodeDistance(i,j,tree);
						dis2 = NodeDistance(i,k,tree);
						dis3 = NodeDistance(i,l,tree);
						dis4 = NodeDistance(j,k,tree);
						dis5 = NodeDistance(j,l,tree);
						dis6 = NodeDistance(k,l,tree);						
						if((dis2+dis5)-(dis1+dis6)>0){
							FindQuartet(x1, x2, x3, x4, location);
							quartet[location[0]][location[1]]++;
						}else if((dis1+dis6)-(dis2+dis5)>0){
							FindQuartet(x1, x3, x2, x4, location);
							quartet[location[0]][location[1]]++;
						}else if((dis1+dis6)-(dis3+dis4)>0){
							FindQuartet(x1, x4, x2, x3, location);
							quartet[location[0]][location[1]]++;
						}else{ /*polytomy*/
							continue;
						}
					}
				}
			}
		}
	}

	return NO_ERROR;
}

int NodeDistance (int node1, int node2, Tree *tree){
	int i, j, nodedistance=1;

	for(i=0; i<tree->numnodes; i++){
		if(tree->postorder[i] == node1){
			for(j=i+1; j<tree->numnodes; j++){
				if(tree->postorder[j] == node2){
					return nodedistance;
				}				
				if(tree->postorder[j] > tree->ntaxa){
					nodedistance++;
				}
			}
		}else if(tree->postorder[i] == node2){
			for(j=i+1; j<tree->numnodes; j++){
				if(tree->postorder[j] == node1){
					return nodedistance;
				}				
				if(tree->postorder[j] > tree->ntaxa-1){/*count # of internal nodes*/
					nodedistance++;
				}
			}
		}else{
			continue;
		}
	}

	return NO_ERROR;
}

int TriplesFreqInatreeCollapse (int **triple, Tree *tree){/*we do NOT need this function. TreeBranchCollapse can collapse short branches*/
	int i, j, k, w, son0, son1, father, cnode;
	int offsprings0[NTAXA], offsprings1[NTAXA], offsprings2[NTAXA], n0, n1, n2;
	int location[2];
	double brlens, minbrlens = collapseBrlens;

	for(i=tree->ntaxa; i<2*tree->ntaxa-1; i++){
		if(i != tree->root){
			brlens = tree->nodes[i].brlens;
			n0 = n1 = 0;
			for(j=0; j<tree->ntaxa; j++){
				offsprings0[j] = 0;
				offsprings1[j] = 0;
			}

			/*find the descendant species for son0 and son1*/ 
			son0 = tree->nodes[i].sons[0];
			son1 = tree->nodes[i].sons[1];
			FindOffsprings (offsprings0, tree, son0);
			FindOffsprings (offsprings1, tree, son1);
			for(j=0; j<tree->ntaxa; j++){
				if (offsprings0[j] == 1) offsprings0[n0++] = tree->nodes[j].namenumber;
				if (offsprings1[j] == 1) offsprings1[n1++] = tree->nodes[j].namenumber;
			}

			cnode = i;
			father = tree->nodes[i].father;
			
			while (father != -1){
				/*if the triples are ignored if brlens < minbrlens, i.e., collapse*/	
				if(brlens > minbrlens){
					for(j=0; j<tree->ntaxa; j++){
						offsprings2[j] = 0;
					} 
							
					n2 = 0;

					if(tree->nodes[father].sons[0] == cnode){
						FindOffsprings (offsprings2, tree, tree->nodes[father].sons[1]);
					}else{
						FindOffsprings (offsprings2, tree, tree->nodes[father].sons[0]);
					}

					for(j=0; j<tree->ntaxa; j++){
						if (offsprings2[j] == 1){
							offsprings2[n2++] = tree->nodes[j].namenumber;
						}
					}
						
					for(j=0; j<n0; j++){
						for(k=0; k<n1; k++){
							/*if(offsprings0[j] == offsprings1[k])
								continue;*/
							for(w=0; w<n2; w++){
								/*if(offsprings0[j] == offsprings2[w] || offsprings1[k] == offsprings2[w])
									continue;*/
								FindTriple (offsprings0[j], offsprings1[k], offsprings2[w], location);
								triple[location[0]][location[1]]++;
							}
						}
					}
				}
				/*updating brlens as it moves up to the father node*/
				brlens += tree->nodes[father].brlens;
				cnode = father;
				father = tree->nodes[father].father;
			}					
		}
	}
	
	return NO_ERROR;
}

void FindOffsprings (int *offsprings, Tree *tree, int inode){
	int i;

	if(inode < tree->ntaxa){
		offsprings[inode] = 1;
	}else{
		for(i=0; i<tree->nodes[inode].nson; i++){
			FindOffsprings (offsprings, tree, tree->nodes[inode].sons[i]);
		}
	}
}

int FindTriple (int n1, int n2, int n3, int *location){
	int i, number=0, node1, node2, node3;

	/*order the node numbers and the topology 0 or 1 or 2 is saved in location[1]*/
	if(n1 < n2 && n2 < n3){
		node1 = n1;
		node2 = n2;
		node3 = n3;
		location[1] = 0;
	}else if(n2 < n1 && n1 < n3){
		node1 = n2;
		node2 = n1;
		node3 = n3;
		location[1] = 0;
	}else if(n1 < n3 && n3 < n2){
		node1 = n1;
		node2 = n3;
		node3 = n2;
		location[1] = 1;
	}else if(n2 < n3 && n3 < n1){
		node1 = n2;
		node2 = n3;
		node3 = n1;
		location[1] = 1;
	}else if(n3 < n1 && n1 < n2){
		node1 = n3;
		node2 = n1;
		node3 = n2;
		location[1] = 2;
	}else if(n3 < n2 && n2 < n1){
		node1 = n3;
		node2 = n2;
		node3 = n1;
		location[1] = 2;
	}

	/*calculate the location of the triple (node1, node2, node3) in triplematrix*/
	for(i=1; i<=node1; i++){
		number += (sptree.ntaxa-i)*(sptree.ntaxa-i-1)/2;
	}
	for(i=node1+2; i<=node2; i++){
		number += (sptree.ntaxa-i);
	}
	number += (node3 - node2 -1);
	location[0] = number;

	return NO_ERROR;
}


int FindQuartet (int n1, int n2, int n3, int n4, int *location){
	int i, number=0, nodes[4];

	/*order the node numbers and the topology 0 or 1 or 2 is saved in location[1]*/
	if(n1 > n2){
		nodes[0] = n2;
		nodes[1] = n1;
	}else{
		nodes[0] = n1;
		nodes[1] = n2;
	}

	if(n3 > n4){
		nodes[2] = n4;
		nodes[3] = n3;
	}else{
		nodes[2] = n3;
		nodes[3] = n4;
	}

	if(nodes[0] < nodes[2]){
		if(nodes[1] < nodes[2]){
			location[1] = 0;
		}else if(nodes[1] > nodes[3]){
			location[1] = 2;
		}else{
			location[1] = 1;
		}
	}else{
		if(nodes[3] < nodes[0]){
			location[1] = 0;
		}else if(nodes[3] > nodes[1]){
			location[1] = 2;
		}else{
			location[1] = 1;
		}
	}

	/*sort node numbers*/
	if(nodes[0] < nodes[2]){
		if(nodes[1] > nodes[2]){
			number = nodes[1];
			nodes[1] = nodes[2];
			nodes[2] = number;
		}
		if(nodes[2] > nodes[3]){
			number = nodes[2];
			nodes[2] = nodes[3];
			nodes[3] = number;
		}
	}else{
		number = nodes[0];
		nodes[0] = nodes[2];
		nodes[2] = number;
		number = nodes[1];
		nodes[1] = nodes[3];
		nodes[3] = number;
		if(nodes[1] > nodes[2]){
			number = nodes[1];
			nodes[1] = nodes[2];
			nodes[2] = number;
		}
		if(nodes[2] > nodes[3]){
			number = nodes[2];
			nodes[2] = nodes[3];
			nodes[3] = number;
		}
	}

	/*calculate the location of the quartet in quartetmatrix*/
	number = 0;
	for(i=1; i<nodes[0]+1; i++){
		number += (sptree.ntaxa-i)*(sptree.ntaxa-i-1)*(sptree.ntaxa-i-2)/6;
	}
	for(i=nodes[0]+2; i<=nodes[1]; i++){
		number += (sptree.ntaxa-i)*(sptree.ntaxa-i-1)/2;
	}
	for(i=nodes[1]+2; i<=nodes[2]; i++){
		number += (sptree.ntaxa-i);
	}
	number += (nodes[3]-nodes[2]-1);

	location[0] = number;

	//printf("%d %d %d %d %d %d\n",nodes[0],nodes[1],nodes[2],nodes[3],number,location[1]);

	return NO_ERROR;
}

int ReadaTree (FILE *fTree, Tree *tree){
/* 
   Both names and numbers for species are accepted.  
   Species names are considered case-sensitive, with trailing blanks ignored.
*/
   	int cnode, cfather=-1, taxa=0;  /* current node and father */
   	int inodeb=0;  /* node number that will have the next branch length */
   	int i, x, level=0, ch=' ', index=0, treestrlength, ntaxa=1, nleft=0, nright=0;
   	char skips[]="\"\'", treestring[MAX_STRING_LENGTH], str[50];
   	int nnode;   
	
	/*move to the first character of a treestring*/
	while(ch != '('){
      	ch = fgetc(fTree);
    }

	/*read treestring and find ntaxa*/	
	while(ch != ';'){		
		treestring[index++] = ch;
		if(ch == ','){
			ntaxa++;
		}
		if(ch == '('){
			nleft++;
		}	
		if(ch == ')'){
			nright++;
		}	
		ch = fgetc(fTree);	
	}

	/*add ending to treestring*/
	treestring[index++] = ';';
	treestrlength = index;

	/*check if # of ( == # of )*/
	if(nleft != nright){
		printf("# of ( != # of )");
		return ERROR;
	}

	/*read tree parameters*/
	tree->ntaxa = ntaxa;
	tree->numnodes = ntaxa + nright;
	nnode = tree->ntaxa; 

   	for(i=0; i<tree->numnodes; i++) {
      		tree->nodes[i].father = -1;
			tree->nodes[i].brlens = NA;
      		tree->nodes[i].nson = 0; 
			tree->nodes[i].sons[0] = -1;
			tree->nodes[i].sons[1] = -1;
   	}
	
   	for (i=0; i<treestrlength; i++) {
		ch = treestring[i];	
		if (!isgraph(ch) || ch==skips[0] || ch==skips[1]){
			continue;
		}else if (ch=='('){
			level++;
			cnode=nnode++;

			if(nnode > 2*(tree->ntaxa)-1){
					printf("check tree: perhaps too many '('s");
					return ERROR;
			}

			if (cfather>=0) {
					tree->nodes[cfather].sons[tree->nodes[cfather].nson++] = cnode;
					tree->nodes[cnode].father=cfather;
			}else{
				tree->root=cnode;
			}          			
			cfather=cnode;
		}else if (ch==')') { 
			level--;  
			inodeb=cfather; 
			cfather=tree->nodes[cfather].father; 
		}else if (ch==':') {
			index = 0;
			while(ch != ',' && ch != ')' && ch != ';'){
				ch = treestring[++i];
				str[index++] = ch;	
			}
			i--;
			str[--index] = '\0';
			sscanf(str, "%lf", &tree->nodes[inodeb].brlens);			
		}else if (ch==',') {

		}else if (ch==';' && level!=0) {
            printf("; in treefile");
            return ERROR;
        }else if (isdigit(ch)){ 
			index = 0;
			while(isdigit(treestring[i])){
				str[index++] = treestring[i];
				i++;
			}			
			inodeb = atoi(str);
			inodeb--;
			tree->nodes[inodeb].father=cfather;
			tree->nodes[cfather].sons[tree->nodes[cfather].nson++]=inodeb;
		}else if (isalpha(ch)){		
			index = 0;
			while(ch != ':' && ch != ',' && ch != ')'){					
				tree->nodes[taxa].taxaname[index++] = ch;
				ch = treestring[++i];
			}
			/*add the null character at the end of the string*/
			tree->nodes[taxa].taxaname[index] = '\0';
			i--;
			tree->nodes[taxa].father = cfather;
			tree->nodes[cfather].sons[tree->nodes[cfather].nson++] = taxa;
			inodeb = taxa;
			taxa++;
		}
   	}

	/*sort sons for each internal node. This can simplify other calculation if sons are sorted*/
	for(i=tree->ntaxa; i<nnode; i++){
		for(level=0; level<tree->nodes[i].nson-1; level++){
			for(index=level+1; index<tree->nodes[i].nson; index++){
				if(tree->nodes[i].sons[level]>tree->nodes[i].sons[index]){
					x = tree->nodes[i].sons[level];
					tree->nodes[i].sons[level] = tree->nodes[i].sons[index];
					tree->nodes[i].sons[index] = x;
				}
			}
		}
	}

	/*check if the tree is rooted*/
	if(tree->nodes[tree->root].nson == 2){
		tree->isrooted = 1;
	}else{
		tree->isrooted = 0;
	}

	/*double check on the tree*/
	if(tree->numnodes != nnode){
		printf("Errors in ReadaTree\n");
		return ERROR;
	};

	if(DEBUG){
		int j;
		for(i=0; i<tree->numnodes; i++) {
			printf("node %d %d %lf", tree->nodes[i].father, tree->nodes[i].nson, tree->nodes[i].brlens);
			for(j=0; j<tree->nodes[i].nson; j++){
				printf(" %d ", tree->nodes[i].sons[j]);
			}
			printf("\n");
		}
		printf("%s %s %d %d\n", tree->nodes[0].taxaname, tree->nodes[1].taxaname, tree->isrooted, tree->root);	
		PrintTree (tree, tree->root, 1, 1, 0, 0, 1);
		printf("%s", printString);		
		return ERROR;
	}

	
	if((!tree->isrooted) && (calculatedistance == 0 || calculatedistance == 1)){
		printf("This gene tree is not a rooted tree!\n");
		return ERROR;
	}

	/*create postorder tree traversal*/
	postorderindex = 0;
	TreeTraversalPostorder (tree->root, tree);

	/*create preorder tree traversal*/
	postorderindex = 0;
	TreeTraversalPreorder (tree->root, tree);
   
   	return NO_ERROR;
}

int ReadaTreeString (char *treestring, Tree *tree, bool taxaorder){
/* 
   Both names and numbers for species are accepted.  
   Species names are considered case-sensitive, with trailing blanks ignored.
*/
   	int cnode, cfather=-1, taxa=0;  /* current node and father */
   	int inodeb=0;  /* node number that will have the next branch length */
   	int i, j, x, level=0, ch=' ', index=0, treestrlength, ntaxa=1, nleft=0, nright=0;
   	char skips[]="\"\'", str[50];
   	int nnode;   
	
	/*read treestring and find ntaxa*/
	treestrlength = strlen(treestring);	
	for(i=0; i<treestrlength; i++){		
		ch = treestring[i];
		if(ch == ','){
			ntaxa++;
		}
		if(ch == '('){
			nleft++;
		}	
		if(ch == ')'){
			nright++;
		}		
	}

	/*check if # of ( == # of )*/
	if(nleft != nright){
		printf("# of ( != # of )");
		return ERROR;
	}

	/*read tree parameters*/
	tree->ntaxa = ntaxa;
	tree->numnodes = ntaxa + nright;
	nnode = tree->ntaxa; 

   	for(i=0; i<tree->numnodes; i++) {
      		tree->nodes[i].father = -1;
			tree->nodes[i].brlens = NA;
      		tree->nodes[i].nson = 0; 
			tree->nodes[i].sons[0] = -1;
			tree->nodes[i].sons[1] = -1;
   	}
	
   	for (i=0; i<treestrlength; i++) {
		ch = treestring[i];	
		if (!isgraph(ch) || ch==skips[0] || ch==skips[1]){
			continue;
		}else if (ch=='('){
			level++;
			cnode=nnode++;

			if(nnode > 2*(tree->ntaxa)-1){
					printf("check tree: perhaps too many '('s");
					return ERROR;
			}

			if (cfather>=0) {
					tree->nodes[cfather].sons[tree->nodes[cfather].nson++] = cnode;
					tree->nodes[cnode].father=cfather;
			}else{
				tree->root=cnode;
			}          			
			cfather=cnode;
		}else if (ch==')') { 
			level--;  
			inodeb=cfather; 
			cfather=tree->nodes[cfather].father; 
		}else if (ch==':') {
			index = 0;
			while(ch != ',' && ch != ')' && ch != ';'){
				ch = treestring[++i];
				str[index++] = ch;	
			}
			i--;
			str[--index] = '\0';
			sscanf(str, "%lf", &tree->nodes[inodeb].brlens);			
		}else if (ch==',') {

		}else if (ch==';' && level!=0) {
            printf("; in treefile");
            return ERROR;
        }else if (isdigit(ch)){ 
			index = 0;
			while(isdigit(treestring[i])){
				str[index++] = treestring[i];
				i++;
			}			
			inodeb = atoi(str);
			inodeb--;
			tree->nodes[inodeb].father=cfather;
			tree->nodes[cfather].sons[tree->nodes[cfather].nson++]=inodeb;
		}else if (isalpha(ch)){		
			index = 0;
			while(ch != ':' && ch != ',' && ch != ')'){					
				str[index++] = ch;
				ch = treestring[++i];
			}
			/*add the null character at the end of the string*/
			str[index] = '\0';
			/*read taxa name*/
			if(taxaorder){
				/*find taxon number*/
				for(j=0; j<ntaxa; j++){
					/*global variable sptable.speciesname*/
					if(!strcmp(str, sptable.speciesname[j])){
						taxa = j;
						break;
					}
				}
				if(j == ntaxa){
					printf("cannot find species %s\n", str);
					return ERROR;
				}
				strcpy(tree->nodes[taxa].taxaname,str);
				i--;
				tree->nodes[taxa].father = cfather;
				tree->nodes[cfather].sons[tree->nodes[cfather].nson++] = taxa;
				inodeb = taxa;
			}else{
				strcpy(tree->nodes[taxa].taxaname,str);
				i--;
				tree->nodes[taxa].father = cfather;
				tree->nodes[cfather].sons[tree->nodes[cfather].nson++] = taxa;
				inodeb = taxa;
				taxa++;
			}
		}
   	}

	/*sort sons for each internal node. This can simplify other calculation if sons are sorted*/
	for(i=tree->ntaxa; i<nnode; i++){
		for(level=0; level<tree->nodes[i].nson-1; level++){
			for(index=level+1; index<tree->nodes[i].nson; index++){
				if(tree->nodes[i].sons[level]>tree->nodes[i].sons[index]){
					x = tree->nodes[i].sons[level];
					tree->nodes[i].sons[level] = tree->nodes[i].sons[index];
					tree->nodes[i].sons[index] = x;
				}
			}
		}
	}

	/*check if the tree is rooted*/
	if(tree->nodes[tree->root].nson == 2){
		tree->isrooted = 1;
	}else{
		tree->isrooted = 0;
	}

	/*double check on the tree*/
	if(tree->numnodes != nnode){
		printf("Errors in ReadaTree\n");
		return ERROR;
	};

	if(DEBUG){
		int j;
		for(i=0; i<tree->numnodes; i++) {
			printf("node %d %d %lf", tree->nodes[i].father, tree->nodes[i].nson, tree->nodes[i].brlens);
			for(j=0; j<tree->nodes[i].nson; j++){
				printf(" %d ", tree->nodes[i].sons[j]);
			}
			printf("\n");
		}
		printf("%s %s %d %d\n", tree->nodes[0].taxaname, tree->nodes[1].taxaname, tree->isrooted, tree->root);	
		PrintTree (tree, tree->root, 1, 1, 0, 0, 1);
		printf("%s", printString);		
		return ERROR;
	}

	
	if((!tree->isrooted) && (calculatedistance == 0 || calculatedistance == 1)){
		printf("This gene tree is not a rooted tree!\n");
		return ERROR;
	}

	/*create postorder tree traversal*/
	postorderindex = 0;
	TreeTraversalPostorder (tree->root, tree);

	/*create preorder tree traversal*/
	postorderindex = 0;
	TreeTraversalPreorder (tree->root, tree);
   
   	return NO_ERROR;
}




/*PrinTree prints trees by taxa numbers; inode is the root of the tree*/
int PrintTree (Tree *tree, int inode, int showName, int showBrlens, int showTheta, int showSupport, int isRooted){

	char	*tempStr;
	int   tempStrSize;

	/* allocate the print string */
	printStringSize = 200;
	printString = (char *)malloc(printStringSize * sizeof(char));
	if (!printString){
		Print ("%s   Problem allocating printString (%d)\n", spacer, printStringSize * sizeof(char));
		return ERROR;
	}
	*printString = '\0';

	tempStrSize = 200;
	tempStr = (char *) malloc(tempStrSize * sizeof(char));
	if (!tempStr){
		Print ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));
		return ERROR;
	}

    /*print "(" and population size*/
	SaveSprintf (&tempStr, &tempStrSize,"(");
	AddToPrintString (tempStr);
					
	WriteTreeToFile (tree, tree->root, showName, showBrlens, showTheta, showSupport, isRooted);

	if(showSupport == YES) 
		SaveSprintf (&tempStr, &tempStrSize,")[&concordance=%.2f]", tree->nodes[tree->root].support);
	else 
		SaveSprintf (&tempStr, &tempStrSize,")");

	AddToPrintString (tempStr);

	if(showTheta == YES) 
		SaveSprintf (&tempStr, &tempStrSize,"[&theta=%.2f];\n",tree->nodes[tree->root].theta);
	else 
		SaveSprintf (&tempStr, &tempStrSize,";\n");

	AddToPrintString (tempStr);
	free (tempStr); 

	return NO_ERROR;					
}

/*PrinPhylipTree prints trees by taxa names; inode is the root of the tree*/
int PrintPhylipTree (Tree *tree, int inode, int showBrlens, int showTheta, int isRooted){
	char	*tempStr;
	int     tempStrSize;
	
	/* allocate the print string */
	printStringSize = 200;
	printString = (char *)malloc((size_t) (printStringSize * sizeof(char)));
	if (!printString){
		Print ("%s   Problem allocating printString (%d)\n", spacer, printStringSize * sizeof(char));
		return ERROR;
	}
	*printString = '\0';
	
	tempStrSize = 200;
	tempStr = (char *) malloc((size_t) (tempStrSize * sizeof(char)));
	if (!tempStr){
		Print ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));
		return ERROR;
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
	
	return NO_ERROR;					
}

void WriteTreeToFile (Tree *tree, int inode, int showName, int showBrlens, int showTheta, int showSupport, int isRooted){
	char	*tempStr;
	int     i, tempStrSize = 200;

	tempStr = (char *) malloc((size_t) (tempStrSize * sizeof(char)));
	if (!tempStr){
		Print ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));
	}
	
	if (tree->nodes[inode].nson == 0){
		if(showName == YES){
			SaveSprintf (&tempStr, &tempStrSize, "%s", tree->nodes[inode].taxaname);
		}else{ 
			SaveSprintf (&tempStr, &tempStrSize, "%d", inode+1);
		}
		AddToPrintString (tempStr);

		if (showBrlens == YES){
			SaveSprintf (&tempStr, &tempStrSize, ":%.2f", tree->nodes[inode].brlens);
			AddToPrintString (tempStr);
		}
	}else{
		if (inode != tree->root){
			SaveSprintf (&tempStr, &tempStrSize, "(");
			AddToPrintString (tempStr);
		}

		for(i=0; i<tree->nodes[inode].nson-1; i++){
			WriteTreeToFile (tree,tree->nodes[inode].sons[i],  showName, showBrlens, showTheta, showSupport, isRooted);
			SaveSprintf (&tempStr, &tempStrSize, ",");
			AddToPrintString (tempStr);
		}
		WriteTreeToFile (tree,tree->nodes[inode].sons[tree->nodes[inode].nson-1], showName, showBrlens, showTheta, showSupport, isRooted);	

		if (inode != tree->root){
			if (tree->nodes[inode].father == tree->root && isRooted == NO){
				if (showBrlens == YES){
					SaveSprintf (&tempStr, &tempStrSize, ",%d:%.2f", tree->nodes[inode].father + 1, tree->nodes[tree->nodes[inode].father].brlens);
					AddToPrintString (tempStr);
	
					if((tree->nodes[tree->nodes[inode].father].theta>0) && showTheta == YES) {
						SaveSprintf (&tempStr, &tempStrSize, "[&theta=%.2f]", tree->nodes[tree->nodes[inode].father].theta);
						AddToPrintString (tempStr);
					}
				}else{
					SaveSprintf (&tempStr, &tempStrSize, ",%d", tree->nodes[inode].father + 1);
					AddToPrintString (tempStr);
				}
			}

			SaveSprintf (&tempStr, &tempStrSize, ")");
			AddToPrintString (tempStr);
		
			if (isRooted == YES) /*tree->nodes[inode].father != tree->root)*/{
				if(tree->nodes[inode].support>0 && showSupport == YES){
					SaveSprintf (&tempStr, &tempStrSize,"[&concordance=%.2f]", tree->nodes[inode].support);
					AddToPrintString (tempStr);
				}
				if(tree->nodes[inode].brlens>0 && showBrlens == YES){
					SaveSprintf (&tempStr, &tempStrSize,":%.2f", tree->nodes[inode].brlens);
					AddToPrintString (tempStr);
				}
				if((tree->nodes[inode].theta > 0) && showTheta == YES){
					SaveSprintf (&tempStr, &tempStrSize, "[&theta=%.2f]", tree->nodes[inode].theta);
					AddToPrintString (tempStr);
				}
			}				
		}
	}
	free (tempStr);	
}

void WritePhylipTreeToFile (Tree *tree, int inode, int showBrlens, int showTheta, int isRooted){	
	char			*tempStr;
	int             tempStrSize = 200;
	
	tempStr = (char *) malloc((size_t) (tempStrSize * sizeof(char)));
	if (!tempStr){
		Print ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));
	}
	
	if (tree->nodes[inode].nson == 0){
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
	}else{
		if (inode != tree->root){
			SaveSprintf (&tempStr, &tempStrSize, "(");
			AddToPrintString (tempStr);
		}
		WritePhylipTreeToFile (tree,tree->nodes[inode].sons[0],  showBrlens, showTheta, isRooted);
		SaveSprintf (&tempStr, &tempStrSize, ",");
		AddToPrintString (tempStr);
		WritePhylipTreeToFile (tree,tree->nodes[inode].sons[1], showBrlens, showTheta, isRooted);	
		if (inode != tree->root){
			if (tree->nodes[inode].father == tree->root && isRooted == NO){
				if (showBrlens == YES){
					SaveSprintf (&tempStr, &tempStrSize, ",%s:%lf", tree->nodes[tree->nodes[inode].father].taxaname, tree->nodes[tree->nodes[inode].father].brlens);
					AddToPrintString (tempStr);
					
					if((tree->nodes[tree->nodes[inode].father].theta>0) && showTheta == YES) {
						SaveSprintf (&tempStr, &tempStrSize, "[#%lf]", tree->nodes[tree->nodes[inode].father].theta);
						AddToPrintString (tempStr);
					}
				}else{
					SaveSprintf (&tempStr, &tempStrSize, ",%s", tree->nodes[tree->nodes[inode].father].taxaname);
					AddToPrintString (tempStr);
				}
			}
			
			if (showBrlens == YES && isRooted == YES) /*tree->nodes[inode].father != tree->root)*/{
				SaveSprintf (&tempStr, &tempStrSize,"):%lf", tree->nodes[inode].brlens);
				AddToPrintString (tempStr);
				if((tree->nodes[inode].theta > 0) && showTheta == YES)
				{
					SaveSprintf (&tempStr, &tempStrSize, "[#%lf]", tree->nodes[inode].theta);
					AddToPrintString (tempStr);
				}
			}else{
				SaveSprintf (&tempStr, &tempStrSize, ")");
				AddToPrintString (tempStr);
			}					
		}
	}
	free (tempStr);
}

int QuartetDistance (char *treefile, int ntrees){
	char **treenames, **quartetstring, a, b;
	int i, j, k, dist=0;
	long nquartets = sptree.ntaxa*(sptree.ntaxa-1)*(sptree.ntaxa-2)*(sptree.ntaxa-3)/24;
	char quartetfile[100], outfile[100];
	FILE *fQuartet, *foutfile;

	/*allocate memory for treenames and quartetstring*/
	treenames = (char**)calloc(ntrees, sizeof(char*));
    treenames[0] = (char*)calloc(50*ntrees, sizeof(char));
   	for(i=0; i<ntrees; i++){
   		treenames[i] = treenames[0] + i*50;
	}
	
	quartetstring = (char**)calloc(ntrees, sizeof(char*));
    quartetstring[0] = (char*)calloc((nquartets+1)*ntrees, sizeof(char));
   	for(i=0; i<ntrees; i++){
   		quartetstring[i] = quartetstring[0] + i*(nquartets+1);
	}

	/*read treenames and quartetstring*/
	sprintf(quartetfile, "%s_quartet.list", treefile);
	sprintf(outfile, "%s_quartet.dis", treefile);
	fQuartet = fopen(quartetfile,"r");
	foutfile = fopen(outfile,"w");

	for(i=0; i<ntrees; i++){
		if(fscanf(fQuartet, "%s%s", treenames[i], quartetstring[i]) != 2){
			printf("Error: cannot read Quartet file %s\n", quartetfile);
			return ERROR;
		}
	}

	/*triple distance*/
	for(i=0; i<ntrees; i++){
		fprintf(foutfile,"%s\t", treenames[i]);
		for(j=0; j<ntrees; j++){
			dist = 0;
			for(k=0; k<nquartets; k++){				
				a = quartetstring[i][k];
				b = quartetstring[j][k];
				if(a>b || a<b){
					dist++;
				}				
			}
			fprintf(foutfile,"%d\t", dist);
		}
		fprintf(foutfile,"\n");
	}

	/*free memory*/	
	free(treenames[0]);
	free(treenames);
	free(quartetstring[0]);
	free(quartetstring);
	fclose(fQuartet);
	fclose(foutfile);

	return NO_ERROR;
}

int QuartetsList (char *treefile, int ntrees){
	int **genequartetmatrix, x1, x2, x3;
	int nquartets = sptree.ntaxa*(sptree.ntaxa-1)*(sptree.ntaxa-2)*(sptree.ntaxa-3)/24;
	int i, j, k, stop;
	Tree genetree;
	char outfile[100];
	FILE *fTree, *foutfile;

	sprintf(outfile, "%s_quartet.list", treefile);
	fTree = fopen(treefile,"r");
	foutfile = fopen(outfile,"w");

	/*allocate memory for genequartetmatrix*/
	genequartetmatrix = (int**)calloc(nquartets, sizeof(int*));
    genequartetmatrix[0] = (int*)calloc(4*nquartets, sizeof(int));
   	for(i = 0; i < nquartets; i++){
		genequartetmatrix[i] = genequartetmatrix[0] + i*4;
	}
  	if(!genequartetmatrix){
		printf("allocating problems for genequartetmatrix!\n");
	   	return ERROR;
	} 

	/*find gene tree quartet*/
	for(i=0; i<ntrees; i++){
		if(ReadaTree(fTree, &genetree) == ERROR) {
			printf("Errors in the gene tree %d; It must be a rooted tree.\n", i+1);
			return ERROR;
		}

		/*find gene tree namenumber*/
		for(j=0; j<genetree.ntaxa; j++){
			stop = 0;
			for(k=0; k<totaltaxa; k++){
				if(!strcmp(genetree.nodes[j].taxaname, sptable.taxaname[k])){
					stop = 1; 
					genetree.nodes[j].namenumber = sptable.taxanodenumber[k];
					break;
				}
			}
			if(stop == 0){
				printf("The taxon %s in gene tree %d is missing in the species-allele table in the control file\n", genetree.nodes[j].taxaname, i);
				return ERROR;
			}
		}

		/*reset genequartetmatrix*/
		for(j=0; j<nquartets; j++){
			for(k=0; k<3; k++){
				genequartetmatrix[j][k] = 0;
			}
		}

		if(QuartetsFreqInatree (genequartetmatrix, &genetree) == ERROR){
			printf("Errors in the QuartetsFreqInatree function\n");
			return ERROR;
		}
			
		fprintf(foutfile,"tree%d\t",i+1);
		for(j=0; j<nquartets; j++){
			x1 = genequartetmatrix[j][0];
			x2 = genequartetmatrix[j][1];
			x3 = genequartetmatrix[j][2];
			if(x1>x2 && x1>x3){
				fprintf(foutfile,"1");
			}else if(x2>x1 && x2>x3){
				fprintf(foutfile,"2");
			}else if(x3>x1 && x3>x2){
				fprintf(foutfile,"3");
			}else{/*missing or polytomy*/
				fprintf(foutfile,"0");
			}	
		}
		fprintf(foutfile,"\n");
	}

	free(genequartetmatrix[0]);
	free(genequartetmatrix);
	fclose(fTree);
	fclose(foutfile);
	
	return NO_ERROR;
}

int TriplesList (char *treefile, int ntrees){
	int **genetriplematrix;
	int ntriples = sptree.ntaxa*(sptree.ntaxa-1)*(sptree.ntaxa-2)/6;
	int i, j, k, stop;
	Tree genetree;
	char outfile[100];
	FILE *fTree, *foutfile;

	sprintf(outfile, "%s_triple.list", treefile);

	fTree = fopen(treefile,"r");
	foutfile = fopen(outfile,"w");

	/*allocate memory for genetriplematrix*/
	genetriplematrix = (int**)calloc(ntriples, sizeof(int*));
    genetriplematrix[0] = (int*)calloc(4*sptree.ntaxa*(sptree.ntaxa-1)*(sptree.ntaxa-2)/6, sizeof(int));
   	for(i = 0; i < ntriples; i++){
		genetriplematrix[i] = genetriplematrix[0] + i*4;
	}
  	if(!genetriplematrix){
		printf("allocating problems for genetriplematrix!\n");
	   	return ERROR;
	} 
	
	/*find gene tree triples*/
	for(i=0; i<ntrees; i++){
		if(ReadaTree(fTree, &genetree) == ERROR) {
			printf("Errors in the gene tree %d; It must be a rooted tree.\n", i+1);
			return ERROR;
		}

		printf("Reading gene tree %d ... calculating triple distance ... Done\n", i);

		/*find gene tree namenumber*/
		for(j=0; j<genetree.ntaxa; j++){
			stop = 0;
			for(k=0; k<totaltaxa; k++){
				if(!strcmp(genetree.nodes[j].taxaname, sptable.taxaname[k])){
					stop = 1; 
					genetree.nodes[j].namenumber = sptable.taxanodenumber[k];
					break;
				}
			}
			if(stop == 0){
				printf("The taxon %s in gene tree %d is missing in the species-allele table in the control file\n", genetree.nodes[j].taxaname, i);
				return ERROR;
			}
		}

		/*reset genetriplematrix*/
		for(j=0; j<ntriples; j++){
			for(k=0; k<3; k++){
				genetriplematrix[j][k] = 0;
			}
		}

		if(TriplesFreqInatree (genetriplematrix, &genetree) == ERROR){
			printf("Errors in the TriplesFreqInatree function\n");
			return ERROR;
		}
			
		fprintf(foutfile,"tree%d\t",i+1);
		for(j=0; j<ntriples; j++){
			if(genetriplematrix[j][0] == 1){
				fprintf(foutfile,"1");
			}else if(genetriplematrix[j][1] == 1){
				fprintf(foutfile,"2");
			}else if(genetriplematrix[j][2] == 1){
				fprintf(foutfile,"3");
			}else{
				fprintf(foutfile,"0");
			}
			//fprintf(foutfile,"\t");	
		}
		fprintf(foutfile,"\n");
	}

	free(genetriplematrix[0]);
	free(genetriplematrix);
	fclose(fTree);
	fclose(foutfile);
	
	return NO_ERROR;
}

int TripleDistance (char *treefile, int ntrees){
	char **treenames, **triplestring, a, b;
	int i, j, k, dist=0;
	long ntriples = sptree.ntaxa*(sptree.ntaxa-1)*(sptree.ntaxa-2)/6;
	char triplefile[100], outfile[100];
	FILE *fTriple, *foutfile;

	/*allocate memory for treenames and triplestring*/
	treenames = (char**)calloc(ntrees, sizeof(char*));
    treenames[0] = (char*)calloc(50*ntrees, sizeof(char));
   	for(i=0; i<ntrees; i++){
   		treenames[i] = treenames[0] + i*50;
	}
	
	triplestring = (char**)calloc(ntrees, sizeof(char*));
    triplestring[0] = (char*)calloc((ntriples+1)*ntrees, sizeof(char));
   	for(i=0; i<ntrees; i++){
   		triplestring[i] = triplestring[0] + i*(ntriples+1);
	}

	/*read treenames and triplestring*/
	sprintf(triplefile, "%s_triple.list", treefile);
	sprintf(outfile, "%s_triple.dis", treefile);
	fTriple = fopen(triplefile,"r");
	foutfile = fopen(outfile,"w");

	for(i=0; i<ntrees; i++){
		if(fscanf(fTriple, "%s%s", treenames[i], triplestring[i]) != 2){
			printf("Error: cannot read Triple file %s\n", triplefile);
			return ERROR;
		}
	}

	/*triple distance*/
	for(i=0; i<ntrees; i++){
		fprintf(foutfile,"%s\t", treenames[i]);
		for(j=0; j<ntrees; j++){
			dist = 0;
			for(k=0; k<ntriples; k++){				
				a = triplestring[i][k];
				b = triplestring[j][k];
				if(a>b || a<b){
					dist++;
				}				
			}
			fprintf(foutfile,"%d\t", dist);
		}
		fprintf(foutfile,"\n");
	}

	/*free memory*/	
	free(treenames[0]);
	free(treenames);
	free(triplestring[0]);
	free(triplestring);
	fclose(fTriple);
	fclose(foutfile);

	return NO_ERROR;
}


int TriplesFreqInaBinaryRootedtree (int **triple, Tree *tree){/*this function is replaced by TripleFreqInatree*/
	int i, j, k, w, son0, son1;
	int offsprings0[NTAXA], offsprings1[NTAXA], offsprings2[NTAXA], n0, n1, n2;
	int location[2];
 	
	for(i=tree->ntaxa; i<2*tree->ntaxa-1; i++){
		if(i != tree->root){
			n0 = n1 = n2 = 0;
			for(j=0; j<tree->ntaxa; j++){
				offsprings0[j] = 0;
				offsprings1[j] = 0;
			}

			son0 = tree->nodes[i].sons[0];
			son1 = tree->nodes[i].sons[1];

			/*find the descendant species for son0 and son1*/
			FindOffsprings (offsprings0, tree, son0);
			FindOffsprings (offsprings1, tree, son1);
			
			/*count the number of species*/
			for(j=0; j<tree->ntaxa; j++){
				if(offsprings0[j] == 0 && offsprings1[j] == 0)
					offsprings2[n2++] = tree->nodes[j].namenumber;
				else if (offsprings0[j] == 1 && offsprings1[j] == 0)
					offsprings0[n0++] = tree->nodes[j].namenumber;
				else if (offsprings0[j] == 0 && offsprings1[j] == 1)
					offsprings1[n1++] = tree->nodes[j].namenumber;
				else
					return ERROR;
			}

			for(j=0; j<n0; j++){
				for(k=0; k<n1; k++){
					/*if(offsprings0[j] == offsprings1[k]){
						return ERROR;
					}*/

					for(w=0; w<n2; w++){
						/*if(offsprings0[j] == offsprings2[w] || offsprings1[k] == offsprings2[w]){
							return ERROR;
						}*/
						FindTriple (offsprings0[j], offsprings1[k], offsprings2[w], location);
						triple[location[0]][location[1]]++;
					}
				}
			}
					
		}
	}

	return NO_ERROR;
}













