/*  
 * Joan A. Parera Portell
 * 
 * 15/06/2019
 * 
 * This program performs bootstrap on H-k stacks.
 * 
 * Arguments:
 * program [output file] [H-k RF 1] [H-k RF 2] ... [H-k RF N]
 * 
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/* Initialize variables for pseudo-random congruential generator. 
 * Borland parameters are used */
#define M 4294967296
#define A 22695477
#define C 1
/* Number of steps in previous stacking (hk_sum.c and hkstacking.c) */
#define STEP 100
/* Number of bootstrap iterations */
#define NITER 100
/* Initialize seed */
long long seed = 1;

/*-----------Pseudo-random number generator function [0, max]---------*/
int rnd(int max)
{
	int out, n;
	double f;
	
	if(max < 1)
	{
		printf("Error, MAX value is too low");
	}
		
	/* Linear mixed congruential generator */
	seed = (A*seed+C)%M;
	/* Scale pseudo-random [0,1] number according to max */
	f = (double)seed/M;
	out = f*max;
	
	/* Returns integer */
	return out;
}

/*-----------------------------Bootstrap------------------------------*/
int main(int argc, char **argv)
{
	int n, m, w, r, z, niter;
	float max, h[STEP*STEP], k[STEP*STEP], s[STEP*STEP][argc-2],
	sum[STEP*STEP], hh, kk, ss, avh, avk, hsum, ksum, covh, covk, covhk, 
	hstd, kstd, corr, mh[NITER], mk[NITER];
	char file[150], outfile[150];
	FILE *fout, *fcov;
  
	if(argc < 3)
	{
		printf("Usage: program [output] [H-k RF 1] ... [H-k RF N]");
		exit(0);	
	}
	
	/* "Warm up" pseudo-random generator */
	for(n=0; n<100; n++)
	{
		seed = (A*seed+C)%M;
	}
	
	/* Begin bootstrap iteration */
	printf("Starting Bootstrap. Iteration number:\n");
	for(z=0; z<NITER; z++)
	{
		/* Iterate through random stack files passed as args */
		printf("%d...", z+1);
		max = -2000.0;
		for(n=2; n<argc; n++)
		{
			/* Call rnd to return random int between 2 to N */
			r = rnd(argc-2)+2;
			/* Open file number r */
			strcpy(file, argv[r]);
			fout = fopen(file, "r");
			m = 0;
			w = 3;
			while(m<STEP*STEP && w==3)
			{
				/* Read file and store values */
				w = fscanf(fout, "%f,%f,%f", &hh, &kk, &ss);
				
				/* Add stack values m to file vector n */
				h[m] = hh;
				k[m] = kk;
				s[m][n-2] = ss;
				m++;
			}
			fclose(fout);
		}
		
		/* Stacking */
		for(m=0; m<STEP*STEP; m++)
		{
			sum[m] = 0.0;
			for(n=0; n<argc-2; n++)
			{
				sum[m] = sum[m] + s[m][n];
			}
			
			/* Find if the stack amplitude is maximum. If true, save H 
			 * and kappa */
			if (sum[m] > max)
			{
				max = sum[m];
				mh[z] = h[m];
				mk[z] = k[m];
			}
		}
		printf("%f, %f\n", mh[z], mk[z]);
	}

	/* Obtain average H (avh) and k (avk) over niter iterations */
	hsum = 0.0;
	ksum = 0.0;
	for(z=0; z<NITER; z++)
	{
		hsum = hsum + mh[z];
		ksum = ksum + mk[z];
	}
	avh = hsum/NITER;
	avk = ksum/NITER;
	printf("\nAverage: %f, %f\n", avh, avk);
	
	/* Compute covariance matrix */
	covh = 0.0;
	covk = 0.0;
	covhk = 0.0;
	for(z=0; z<NITER; z++)
	{
		covh = covh+((mh[z]-avh)*(mh[z]-avh));
		covk = covk+((mk[z]-avk)*(mk[z]-avk));
		covhk = covhk+((mh[z]-avh)*(mk[z]-avk));
	}
	
	covh = (float)covh/(NITER-1);
	covk = (float)covk/(NITER-1);
	covhk = (float)covhk/(NITER-1);
	
	/* Compute standard error */
	hstd = sqrtf(covh);
	kstd = sqrtf(covk);
	
	/* Compute correlation */
	corr = covhk/(kstd*hstd);
	
	/* Write results to file */
	strcpy(outfile, argv[1]);
	fcov = fopen(outfile, "w");
	fprintf(fcov, "Hstd,Kstd,CORR\n%f,%f,%f\n", hstd, kstd, corr);
	printf("Hstd, Kstd, CORR: %f, %f, %f\n", hstd, kstd, corr);
	fclose(fcov);
	
	return 0;
}

