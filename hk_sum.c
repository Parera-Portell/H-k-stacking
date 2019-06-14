/*  
 * Joan A. Parera Portell
 * 
 * 14/06/2019
 * 
 * This program stacks and normalizes H-kappa RF.
 * 
 * Arguments:
 * program [output file] [H-k RF 1] [H-k RF 2] ... [H-k RF N]
 * 
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/* 
 * Define number of iteration steps (must be the same than that in 
 * hkstacking.c.) 
 */
#define STEP 101

int main(int argc, char **argv)
{
	int n, m, w;
	float min, max, h[STEP*STEP], k[STEP*STEP], s[STEP*STEP][argc-2],
	sum[STEP*STEP], hh, kk, ss, si;
	char file[150], outfile[150];
	FILE *fout, *fsum;
  
	if(argc < 3)
	{
		printf("Usage: program [output txt] [H-k RF 1] [H-k RF 2] ... [H-k RF N]");
		exit(0);	
	}
	
	min = 20000.0;
	max = -20000.0;
	strcpy(outfile, argv[1]);
	
	/* Iteration through stack files passed as args */
	for(n=2; n<argc; n++)
	{
		strcpy(file, argv[n]);
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
	}
	
	/* Stacking */
	for(m=0; m<STEP*STEP; m++)
	{
		sum[m] = 0;
		for(n=0; n<argc-2; n++)
		{
			sum[m] = sum[m] + s[m][n];
		}
		
		/* Find if the stack value is maximum */
		if (sum[m] > max)
		{
			max = sum[m];
		}
		
		/* Find if the stack value is minimum */
		if (sum[m] < min)
		{
			min = sum[m];
		}
	}
	
	fclose(fout);
	
	/* Normalize stack and write to output file */
	fsum = fopen(outfile, "w");
	for(n=0; n<STEP*STEP; n++)
	{
		si = (sum[n]+sqrtf(powf(min,2.0)))/(max+sqrtf(powf(min,2.0)));
		fprintf(fsum, "%2.2f,%1.3f,%1.6f\n", h[n], k[n], si);
	}
	
	fclose(fsum);

	return 0;
}
