/*  
 * Joan A. Parera Portell
 * 
 * 15/06/2019
 * 
 * This program performs bootstrap on H-k stacks.
 * 
 * Arguments:
 * program [out error] [out iter] [H-k RF 1] [H-k RF 2] ... [H-k RF N]
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
/* Number of bootstrap iterations */
#define NITER 100
/* Maximum array size */
#define MAX 50000
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
		exit(0);
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
	int n, m, w, r, z, niter, size;
	float max, h[MAX], k[MAX], s[MAX], sum[MAX], hh, kk, ss, avh, avk, 
	hsum, ksum, varh, vark, covhk, hstd, kstd, corr, mh[NITER], mk[NITER];
	char file[150], outfile[150], outiter[150];
	FILE *fout, *fcov, *fiter;
  
	if(argc < 4)
	{
		printf("Usage: program [out error] [out iter] [H-k RF 1] ...");
		exit(0);	
	}
	
	/* "Warm up" pseudo-random generator */
	for(n=0; n<100; n++)
	{
		seed = (A*seed+C)%M;
	}
	
	/* Begin bootstrap iteration */
	printf("\nStarting Bootstrap. Iteration number:\n");
	for(z=0; z<NITER; z++)
	{
		/* Iterate through random stack files passed as args */
		printf("%d,", z+1);
		max = -2000.0;
		for(n=3; n<argc; n++)
		{
			/* Call rnd to return random integer identifier */
			r = rnd(argc-3)+3;
			/* Open file number r */
			strcpy(file, argv[r]);
			fout = fopen(file, "r");
			size = 0;
			w = 3;
			while(w==3)
			{
				/* Read file and store values */
				w = fscanf(fout, "%f,%f,%f", &hh, &kk, &ss);
				
				if (n == 3)
				{
					s[size] = 0;
				}
				
				/* Stacking */
				h[size] = hh;
				k[size] = kk;
				s[size] += ss;
				
				/* Find if the stack value is maximum */
				if (s[size] > max)
				{
					max = s[size];
					mh[z] = h[size];
					mk[z] = k[size];
				}

				size++;
			}
			
			fclose(fout);
		}

		printf("%f,%f\n", mh[z], mk[z]);
	}

	/* Obtain average H (avh) and k (avk) over niter iterations */
	hsum = 0.0;
	ksum = 0.0;
	strcpy(outiter, argv[2]);
	fiter = fopen(outiter, "w");
	for(z=0; z<NITER; z++)
	{
		hsum = hsum + mh[z];
		ksum = ksum + mk[z];
		fprintf(fiter, "%f,%f\n", mh[z], mk[z]);
	}
	avh = hsum/NITER;
	avk = ksum/NITER;
	printf("\nAverage: %f, %f\n", avh, avk);
	
	/* Compute covariance matrix */
	varh = 0.0;
	vark = 0.0;
	covhk = 0.0;
	for(z=0; z<NITER; z++)
	{
		varh = varh+((mh[z]-avh)*(mh[z]-avh));
		vark = vark+((mk[z]-avk)*(mk[z]-avk));
		covhk = covhk+((mh[z]-avh)*(mk[z]-avk));
	}
	
	varh = (float)varh/(NITER-1);
	vark = (float)vark/(NITER-1);
	covhk = (float)covhk/(NITER-1);
	
	/* Compute standard error */
	hstd = sqrtf(varh);
	kstd = sqrtf(vark);
	
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
