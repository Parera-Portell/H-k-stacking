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

/* Maximum array size */
#define MAX 50000

int main(int argc, char **argv)
{
	int n, m, w, size;
	float min, max, h[MAX], k[MAX], s[MAX], hh, kk, ss, si, esth, estk;
	char file[150], outfile[150];
	FILE *fout, *fsum;
  
	if(argc < 3)
	{
		printf("Usage: program [output] [H-k RF 1] ... [H-k RF N]");
		exit(0);	
	}
	
	min = 20000.0;
	max = -20000.0;
	strcpy(outfile, argv[1]);
	
	/* Iteration through stack files passed as args */
	printf("\nReading H-k RF number:\n");
	for(n=2; n<argc; n++)
	{
		strcpy(file, argv[n]);
		fout = fopen(file, "r");
		printf("%d...", n-1);
		size = 0;
		w = 3;
		while(w==3)
		{
			/* Read file and store values */
			w = fscanf(fout, "%f,%f,%f", &hh, &kk, &ss);
			
			if (n == 2)
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
			}
			
			/* Find if the stack value is minimum */
			if (s[size] < min)
			{
				min = s[size];
			}
			
			size++;
		}
	}
	
	fclose(fout);
	size--;
	
	/* Normalize stack and write to output file */
	fsum = fopen(outfile, "w");
	for(n=0; n<size; n++)
	{
		si = (s[n]+sqrtf(powf(min,2.0)))/(max+sqrtf(powf(min,2.0)));
		fprintf(fsum, "%2.2f,%1.3f,%1.5f\n", h[n], k[n], si);
	}
	printf("\nH-k stacking done!\n");
	fclose(fsum);

	return 0;
}
