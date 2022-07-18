double GetTime(int inode,double t[])
{
  double time=0;
  int k,kson;

  FOR(k,sptree.nodes[inode].nson) 
    {kson=sptree.nodes[inode].sons[k];
     time+=(sptree.nodes[kson].divtime+GetTime(kson,t))/sptree.nodes[inode].nson;
     }

  if(time>0){
  t[pigindex]=time;
  pigindex++;}
  return(time);
}

double LnP1 (double t, double l, double m, double r)

{

	double		p0t;
	
	p0t = r*(l-m) / (r*l + (l*(1.0-r)-m)*exp((m-l)*t) );
	
	return (log(1.0/r) + 2.0*log(p0t) + (m-l)*t);

}





double LnVt (double t, double l, double m, double r)

{

	double		p0t;
	
	p0t = r*(l-m) / (r*l + (l*(1.0-r)-m)*exp((m-l)*t) );
	
	return (log(1.0 - (1.0/r) * p0t * exp((m-l)*t)));

}

int LnBirthDeathPriorPr (double *prob, double sR, double eR, double sF)

{

	int			i, j;	
        double                  nt[NS],rootTime;
	

	/* rescale all of the node times on the tree */
		pigindex=0;
        GetTime(sptree.root,nt);rootTime=nt[com.nspecies-2];

	for (i=0; i<com.nspecies-1; i++)
		nt[i] /= rootTime;
		
	/* I think this is correct. It looks as if Yang and Rannala (1997)
	   have the root time constrained to be 1.0. */
	rootTime = 1.0;
							
	/* calculate probabilities of tree */
	if (sR != eR)
		{
		(*prob) = (com.nspecies - 2.0) * log(sR);
		for (i=0; i<com.nspecies-2; i++)
			(*prob) += LnP1 (nt[i], sR, eR, sF) - LnVt (rootTime, sR, eR, sF);
		}
	else
		{
		(*prob) = 0.0;
		for (i=0; i<com.nspecies-2; i++)
			(*prob) += log (1.0 + sF * eR) - (2.0 * log(1.0 + sF * eR * nt[i]));
		}
		
	
	return (NO_ERROR);
		
}



double LnBirthDeathPr_diff(double sR, double eR, double sF,double newdistance,double olddistance,double oldroottime,double newroottime,int change)
{
 if(sR!=eR){
   if(oldroottime==newroottime) return(LnP1(newdistance,sR,eR,sF)-LnP1(olddistance,sR,eR,sF));
   else{
         if(change==0) return(com.treeHeightExp*(newroottime-oldroottime)+(com.nspecies-2)*(LnVt (oldroottime, sR, eR, sF)-LnVt (newroottime, sR, eR, sF)));
         else if(change==1) return(com.treeHeightExp*(newroottime-oldroottime)+(com.nspecies-2)*(LnVt (oldroottime, sR, eR, sF)-LnVt (newroottime, sR, eR, sF))+(LnP1(newdistance,sR,eR,sF)-LnP1(newroottime,sR,eR,sF)));
         else if(change==2) return(com.treeHeightExp*(newroottime-oldroottime)+(com.nspecies-2)*(LnVt (oldroottime, sR, eR, sF)-LnVt (newroottime, sR, eR, sF))+(LnP1(oldroottime,sR,eR,sF)-LnP1(olddistance,sR,eR,sF)));
         else  { printf("there is no this option"); exit(1);}
       }}

 else{ if(newroottime==oldroottime) return(2.0 * log(1.0 + sF * eR * olddistance)-2.0 * log(1.0 + sF * eR * newdistance));
       else{
             if(change==0) return(com.treeHeightExp*(newroottime-oldroottime));
             else if(change==1) return(com.treeHeightExp*(newroottime-oldroottime)+2.0 * log(1.0 + sF * eR * newroottime)-2.0 * log(1.0 + sF * eR * newdistance));
             else if(change==2) return(com.treeHeightExp*(newroottime-oldroottime)+2.0 * log(1.0 + sF * eR * olddistance)-2.0 * log(1.0 + sF * eR * oldroottime));
             else  { printf("there is no this option"); exit(1);}
           }
     }
}


double Move_Extinction (void)

{

	/* change extinction rate using sliding window */
	
	int			isMPriorExp, isValidM;
	double		oldM, newM, window, minM, maxM,ran, sR, eR, sF, oldLnPrior, newLnPrior,lnProposalRatio,lnPriorRatio;


	/* get size of window, centered on current mu value */
	window = mcmc.mvp_eR;

	
	
	/* get minimum and maximum values for mu */
	
		minM = mcmc.extinction[0];
		if (mcmc.sR > mcmc.extinction[1])
			maxM = mcmc.extinction[1];
		else
			maxM = mcmc.sR;
		if (maxM < minM)
			minM = 0.0;
		
	

	/* get old value of mu */
	newM = oldM = mcmc.eR;

	/* change value for mu */
	ran = rndu();
	newM = oldM + window * (ran - 0.5);
	
	/* check that new value is valid */
	isValidM = NO;
	do
		{
		if (newM < minM)
			newM = 2* minM - newM;
		else if (newM > maxM)
			newM = 2 * maxM - newM;
		else
			isValidM = YES;
		} while (isValidM == NO);

	/* get proposal ratio */
	lnProposalRatio = 0.0;
	
	/* calculate prior ratio */
	
	sF = mcmc.sF;
	eR = oldM;
	if (LnBirthDeathPriorPr (&oldLnPrior, sR, eR, sF) == ERROR)
		{
		printf ("   Problem calculating prior for birth-death process\n");
		
		}
	eR = newM;
	if (LnBirthDeathPriorPr (&newLnPrior, sR, eR, sF) == ERROR)
		{
		printf ("  Problem calculating prior for birth-death process\n");
		
		}
	
	lnPriorRatio = newLnPrior - oldLnPrior;
	
	if(lnPriorRatio>0||rndu()<exp(lnPriorRatio)) {mcmc.eR=newM;return(lnPriorRatio);}
                                
	else return (0);

}


double Move_Speciation (void)

{

	/* change speciation rate using sliding window */
	
	int			isLPriorExp, isValidL;
	double		oldL, newL, window, minL, maxL, lambdaExp, ran, sR, eR, sF, oldLnPrior, newLnPrior,lnPriorRatio;
	double		lnProposalRatio;

	/* get size of window, centered on current lambda value */
	window = mcmc.mvp_sR;

	
	/* get minimum and maximum values for lambda */
	
		
		minL = mcmc.speciation[0];
		maxL = mcmc.speciation[1];
	
	
	/* get old value of lambda */
	newL = oldL = mcmc.sR;

	/* change value for lambda */
	ran = rndu();
	newL = oldL + window * (ran - 0.5);
	
	/* check that new value is valid */
	isValidL = NO;
	do
		{
		if (newL < minL)
			newL = 2* minL - newL;
		else if (newL > maxL)
			newL = 2 * maxL - newL;
		else
			isValidL = YES;
		} while (isValidL == NO);

	/* get proposal ratio */
	lnProposalRatio = 0.0;
	
	/* calculate prior ratio */
	
	sF = mcmc.sF;
	sR = oldL;
	if (LnBirthDeathPriorPr (&oldLnPrior, sR, eR, sF) == ERROR)
		{
		printf ("  Problem calculating prior for birth-death process\n");
		
		}
	sR = newL;
	if (LnBirthDeathPriorPr ( &newLnPrior, sR, eR, sF) == ERROR)
		{
		printf ("  Problem calculating prior for birth-death process\n");
		
		}
	
	lnPriorRatio = newLnPrior - oldLnPrior;
	
	
	if(lnPriorRatio>0||rndu()<exp(lnPriorRatio)) {mcmc.sR=newL;return(lnPriorRatio);}
                                
	else return (0);
	
}

