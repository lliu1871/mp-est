#include "mpest.h"

FILE *gfopen(char *filename, char *mode)
{
   FILE *fp=(FILE*)fopen(filename, mode);
   if(fp==NULL) {
      printf("\nerror when opening file %s\n", filename);
      exit(-1);
   }
   return(fp);
}

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

int AddToPrintString (char *tempStr){
	size_t len1, len2;

	len1 = (int) strlen(printString);
	len2 = (int) strlen(tempStr);
	if (len1 + len2 + 5 > printStringSize){
		printStringSize += len1 + len2 - printStringSize + 200;
		printString = realloc((void *)printString, printStringSize * sizeof(char));
		if (!printString){
			Print ("Problem reallocating printString (%d)\n", printStringSize * sizeof(char));
			return ERROR;
		}
	}
	strcat(printString, tempStr);

	return NO_ERROR;
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
		return ERROR;
	}

	rewind(fTree);
	return NO_ERROR;
}

/*///////////////////////////////////////////*/
/*      math function                       */
/*/////////////////////////////////////////*/
long factorial(int n){
   long f, i;
   if (n>10) printf("n>10 in factorial");
   for (i=2,f=1; i<=(long)n; i++) f*=i;
   return (f);
}


void SetSeed (unsigned int seed){
   srand(seed);
}

double rndu (void){
	return (double)rand() / (double)((unsigned)RAND_MAX + 1);
}

double rndgamma1 (double s);
double rndgamma2 (double s);


double rndgamma (double s){
/* random standard gamma (Mean=Var=s,  with shape parameter=s, scale para=1)
      r^(s-1)*exp(-r)
   J. Dagpunar (1988) Principles of random variate generation,
   Clarendon Press, Oxford
   calling rndgamma1() if s<1 or
           rndgamma2() if s>1 or
           exponential if s=1
*/
   double r=0;

   if (s<=0)      puts ("jgl gamma..");
   else if (s<1)  r=rndgamma1 (s);
   else if (s>1)  r=rndgamma2 (s);
   else           r=-log(rndu());
   return (r);
}




double rndgamma1 (double s)
{
/* random standard gamma for s<1
   switching method
*/
   double r, x=0,small=1e-37,w;
   static double a,p,uf,ss=10,d;


   if (s!=ss) {
      a=1-s;
      p=a/(a+s*exp(-a));
      uf=p*pow(small/a,s);
      d=a*log(a);
      ss=s;
   }
   for (;;) {
      r=rndu();
      if (r>p)        x=a-log((1-r)/(1-p)), w=a*log(x)-d;
      else if (r>uf)  x=a*pow(r/p,1/s), w=x;
      else            return (0);
      r=rndu ();
      if (1-r<=w && r>0)
         if (r*(w+1)>=1 || -log(r)<=w)  continue;
      break;
   }
   return (x);
}


double rndgamma2 (double s)
{
/* random standard gamma for s>1
   Best's (1978) t distribution method
*/
   double r,d,f,g,x;
   static double b,h,ss=0;
   if (s!=ss) {
      b=s-1;
      h=sqrt(3*s-0.75);
      ss=s;
   }
   for (;;) {
      r=rndu ();
      g=r-r*r;
      f=(r-0.5)*h/sqrt(g);
      x=b+f;
      if (x <= 0) continue;
      r=rndu();
      d=64*r*r*g*g*g;
      if (d*x < x-2*f*f || log(d) < 2*(b*log(x/b)-f))  break;
   }
   return (x);
}

double IncompleteGamma (double x, double alpha, double ln_gamma_alpha)
{
/* returns the incomplete gamma ratio I(x,alpha) where x is the upper 
           limit of the integration and alpha is the shape parameter.
   returns (-1) if in error
   ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
   (1) series expansion     if (alpha>x || x<=1)
   (2) continued fraction   otherwise
   RATNEST FORTRAN by
   Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
   19: 285-287 (AS32)
*/
   int i;
   double p=alpha, g=ln_gamma_alpha;
   /* double accurate=1e-8, overflow=1e30; */
   double accurate=1e-10, overflow=1e60;
   double factor, gin=0, rn=0, a=0,b=0,an=0,dif=0, term=0, pn[6];


   if (x==0) return (0);
   if (x<0 || p<=0) return (-1);


   factor=exp(p*log(x)-x-g);   
   if (x>1 && x>=p) goto l30;
   /* (1) series expansion */
   gin=1;  term=1;  rn=p;
 l20:
   rn++;
   term*=x/rn;   gin+=term;


   if (term > accurate) goto l20;
   gin*=factor/p;
   goto l50;
 l30:
   /* (2) continued fraction */
   a=1-p;   b=a+x+1;  term=0;
   pn[0]=1;  pn[1]=x;  pn[2]=x+1;  pn[3]=x*b;
   gin=pn[2]/pn[3];
 l32:
   a++;  b+=2;  term++;   an=a*term;
   for (i=0; i<2; i++) pn[i+4]=b*pn[i+2]-an*pn[i];
   if (pn[5] == 0) goto l35;
   rn=pn[4]/pn[5];   dif=fabs(gin-rn);
   if (dif>accurate) goto l34;
   if (dif<=accurate*rn) goto l42;
 l34:
   gin=rn;
 l35:
   for (i=0; i<4; i++) pn[i]=pn[i+2];
   if (fabs(pn[4]) < overflow) goto l32;
   for (i=0; i<4; i++) pn[i]/=overflow;
   goto l32;
 l42:
   gin=1-factor*gin;


 l50:
   return (gin);
}


double LnGamma (double x)
{
/* returns ln(gamma(x)) for x>0, accurate to 10 decimal places.
   Stirling's formula is used for the central polynomial part of the procedure.


   Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
   Communications of the Association for Computing Machinery, 9:684
*/
   double f=0, fneg=0, z, lng;
   int nx=(int)x-1;


   if((double)nx==x && nx>0 && nx<10)
      lng=log((double)factorial(nx));
   else {
      if(x<=0) {
         printf("lnGamma not implemented for x<0");
         if((int)x-x==0) { puts("lnGamma undefined"); return(-1); }
         for (fneg=1; x<0; x++) fneg/=x;
         if(fneg<0) printf("strange!! check lngamma");
         fneg=log(fneg);
      }
      if (x<7) {
         f=1;  z=x-1;
         while (++z<7)  f*=z;
         x=z;   f=-log(f);
      }
      z = 1/(x*x);
      lng = fneg+ f + (x-0.5)*log(x) - x + .918938533204673 
             + (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z
                  +.083333333333333)/x;
   }
   return  lng;
}






