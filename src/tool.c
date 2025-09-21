#include "tool.h"

int FindNtaxa (FILE *fTree, int *ntaxa, int ngene);

char	*printString;                /* string for printing to a file                */
size_t	printStringSize;             /* length of printString                        */

static time_t time_start;
void starttime (void){  
   time_start=time(NULL);
}

void strcase (char *str, int direction){
/* direction = 0: to lower; 1: to upper */
   char *p=str;
   if(direction)  while(*p) { *p=(char)toupper(*p); p++; }
   else           while(*p) { *p=(char)tolower(*p); p++; }
}

void Print (char *format, ...){
	va_list  ptr;
	va_start (ptr, format);
	vprintf (format, ptr);
	fflush (stdout);
	va_end (ptr);
}

int AddToPrintString (char *tempStr){
	size_t len1, len2;

	len1 = (int) strlen(printString);
	len2 = (int) strlen(tempStr);
	if (len1 + len2 + 5 > printStringSize){
		printStringSize += len1 + len2 - printStringSize + 200;
		printString = realloc((void *)printString, printStringSize * sizeof(char));
		if (!printString){
			Print ("Problem reallocating printString (%d)\n", printStringSize * sizeof(char));
			return 1;
		}
	}
	strcat(printString, tempStr);

	return 0;
}


#define TARGETLENDELTA (100)

int SaveSprintf(char **target, int *targetLen, char *fmt, ...) {
  va_list    argp;
  int        len, retval;

  va_start(argp, fmt);

#ifdef VISUAL
  len = _vsnprintf(NULL, 0, fmt, argp);
#else
  len = vsnprintf(NULL, 0, fmt, argp);
#endif

  va_end(argp);

  if(len > *targetLen)
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


/*Find the numbers of taxa for all trees in a tree file*/
int FindNtaxa(FILE *fTree, int *ntaxa, int ngene){ 
	char ch;
	int ncomma, i=0;

	ch = fgetc (fTree);
	while (ch != EOF) {
      	ch = fgetc (fTree);
		ncomma = 0;
		while (ch != ';' && ch != EOF){
			ch = fgetc (fTree);
			if(ch == ','){
				ncomma++;
			}		
		}

      	if(ch == ';'){
			ntaxa[i] = ncomma + 1;
        	i++;
		}

		if(i == ngene){
			break;
		}
	}

	if(i < ngene){
		printf("Check on the number of gene trees in the control file. There are only %d gene trees in the tree file!\n", i);
		return 1;
	}

	rewind(fTree);
	return 0;
}

void SetSeed (unsigned int seed){
   srand(seed);
}

double rndu (void){
	return (double)rand() / (double)((unsigned)RAND_MAX + 1);
}
