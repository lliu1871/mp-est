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

/*global variables*/
extern char		*printString;                /* string for printing to a file                */
extern size_t	printStringSize;             /* length of printString                        */

void SetSeed (unsigned int seed);
double rndu (void);
void Print (char *format, ...);
int SaveSprintf(char **target, int *targetLen, char *fmt, ...);
int AddToPrintString (char *tempStr);

