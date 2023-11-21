#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include "bnldev.c"
#include "gammln.c"
#include "ran1.c"

void main(int argc, char **argv)

{

  double mua=1e-7, mum=1e-6, F=100, Fmua, muapmum;
  int ngrowths=4;
  int pinit[4] = { 0, 20,20,20};// Ninit is 2^pinit for each growth except first
  int Ninit = 1;  // Ninit for very first growth
  long long Nmax = (long long)pow(2.0,31);// max population size for every growth
  int nreps = 125;   // number of replications
  long long NAMB[4], Ntotal;
  double cumsum[4];  // cumulative frequency of each type in the population
  FILE *fpout;
  char filename[100];
  int i, igrowth, irep;
  float bottle, r, r2, r3; // bottleneck sampling fraction, random numbers
  float poidev(float,long *), ran1(long *), gamdev(int, long *), bnldev(float, int, long *);
  long seed=-8;

  // set the random number generator seed according to the command line
  // argument, or else it stays at -1
  if (argc>1) seed = -((long)(atof(argv[1])));

  sprintf(filename, "data/luria_%dtrans_%1.2emua_%1.2emum_%dF_%ldseed.out",ngrowths,mua,mum,(int)F,-seed);
  fpout = fopen(filename,"w");
  if (fpout==NULL) { fprintf(stdout,"Unable to open output file\n"); exit(1);}

  Fmua= F*mua;
  muapmum = mua+mum;
  cumsum[3] = 1;
  for (irep=0;irep<nreps;irep++) {

    NAMB[0] = Ninit;
    NAMB[1] = 0; NAMB[2] = 0; NAMB[3] = 0;
    for (igrowth=0;igrowth<ngrowths;igrowth++) {
      if (igrowth>0) {   // bottleneck before growing except for first growth phase
	  bottle = (float)(pow(2.0,pinit[igrowth]))/Ntotal;
          for (i=0;i<4;i++)
	    NAMB[i] = (long long)bnldev(bottle, NAMB[i], &seed);
	}
        Ntotal = NAMB[0]+NAMB[1]+NAMB[2]+NAMB[3];
    if (Ntotal<=0){
        fprintf(stdout,"Error, total popn size negative on rep %d\n",irep);
        Ntotal=Nmax;
        igrowth = ngrowths;
        }
    else{
	cumsum[0] = (double)NAMB[0]/Ntotal;
	for (i=1;i<3;i++) cumsum[i] = cumsum[i-1] + (double)NAMB[i]/Ntotal;
        while (Ntotal<Nmax) {   // loop to accomplish the growth phase, only daughters mutate
          r = ran1(&seed);
	  if (r<cumsum[0])  {  // N reproduces
            NAMB[0]--;
	    if ((r2 = r/cumsum[0]) > muapmum) NAMB[0]++; // check this first b/c it happens most often
	    else if (r2 > mua) NAMB[2]++;
	    else NAMB[1]++;
        if ((r3 = ran1(&seed)) > muapmum) NAMB[0]++; // check this first b/c it happens most often
	    else if (r3 > mua) NAMB[2]++;
	    else NAMB[1]++;
	  }
	  else if (r<cumsum[1]) { // A reproduces
            NAMB[1]--;
	    if ((r2=ran1(&seed)) > 2*mum*(1-mum)) NAMB[1]=NAMB[1]+2;
	    else if (r2> mum*mum) {NAMB[1]++; NAMB[3]++;}
	    else NAMB[3]=NAMB[3]+2;
	  }
	  else if (r<cumsum[2]) { // M reproduces
            NAMB[2]--;
	    if (r2=ran1(&seed) > 2*Fmua*(1-Fmua)) NAMB[2]=NAMB[2]+2;
	    else if (r2> Fmua*Fmua) {NAMB[2]++; NAMB[3]++;}
	    else NAMB[3]=NAMB[3]+2;
	  }
	  else { // B reproduces
	    NAMB[3]++;
	  }
          Ntotal = Ntotal + 1;
	  cumsum[0] = (double)NAMB[0]/Ntotal;
	  for (i=1;i<3;i++) cumsum[i] = cumsum[i-1] + (double)NAMB[i]/Ntotal;
        } // end of while loop, population has reached Nmax
    } // end of else, population size is still positive
	fprintf(fpout," %lld %lld %lld %lld ",NAMB[0],NAMB[1],NAMB[2],NAMB[3]);
	//fprintf(stdout,"Error, Nmax %lld\n",Nmax);
    } // end of loop on growths, all the growth phases are finished
	fprintf(fpout,"\n");
	fprintf(stdout," %d",irep);
  }  // end of loop on reps
  fclose(fpout);
}

