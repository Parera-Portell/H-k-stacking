/*  
 * Joan A. Parera Portell
 * 
 * 13/06/2019
 * 
 * This program converts a time-domain receiver function (RF) to H-kappa 
 * domain. Use hk_sum later to perform stacking. 
 * 
 * Arguments:
 * program [RF file] [output file] [minimum H] [maximum H] [minimum
 * kappa] [maximum kappa] [RF sample rate] [P wave velocity] [weight 1]
 * [weight 2] [weight 3]
 * 
 * Additional parameters can be set below as constants.
 * 
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sacio.h>

/* Maximum length of the data array */
#define MAX 2048
/* Seconds before time 0 in RF */
#define SEC 10
/* Number of iterations */
#define STEP 101
/* Variable containing ray parameter in SAC header */
#define RAYP "USER8"

int main(int argc, char **argv)
{
	float array[MAX], beg, del, p, hmin, hmax, kmin, kmax, vp, vs, 
	kappa, tps, tppps, tppss, w1, w2, w3, si, sn, vpterm,
	vsterm, deltah, deltak;
	int freq, nlen, nerr, mx = MAX, start, n, m, dtps, dtppps, dtppss;
	char file[150], output[150];
	FILE *fout;
	/* Define 3D array that will contain the H-k RF */
	float h[STEP], k[STEP], s[STEP][STEP];
  
	if(argc != 12)
	{
		printf("Usage: hkstacking file file_out hmin hmax kmin kmax freq vp w1 w2 w3");
		exit(0);	
	}

	/* Copy the name of the file to be read into variable */
	strcpy(file , argv[1]);
	/* Output file */
	strcpy(output , argv[2]);
	/* Minimum H */
	hmin = atof(argv[3]);
	/* Maximum H */
	hmax = atof(argv[4]);
	/* Minimum kappa */
	kmin = atof(argv[5]);
	/* Maximum kappa */
	kmax = atof(argv[6]);
	/* Sampling rate (Hz) */
	freq = atoi(argv[7]);
	/* P wave velocity */
	vp = atof(argv[8]);
	/* Weight 1 */
	w1 = atof(argv[9]);
	/* Weight 2 */
	w2 = atof(argv[10]);
	/* Weight 3 */
	w3 = atof(argv[11]);
	/* Define H and kappa increase */
	deltah = (hmax-hmin)/(STEP-1);
	deltak = (kmax-kmin)/(STEP-1);
	/* Define starting point of array */
	start = SEC*freq;
  
	/* Call rsac1 (SAC library) to read sac file. Returns the array
	 * variable. nlen: array length; beg: beggining time; del: delta
	 * or time sampling; mx: MAX; nerr: error return flag; strlen(file):
	 * length of file path.
	 */
	rsac1(file, array, &nlen, &beg, &del, &mx, &nerr, strlen(file));

	/* Check the error status (0=success) */
	if (nerr != 0) 
	{
		printf("Error reading in SAC file: %s\n", file);
		exit (nerr);
	}
	
	/* Call getfhv (SAC library) to get ray parameter p from header */
	getfhv(RAYP, &p , &nerr , strlen(RAYP));
	
	/* Check the error status (0=success) */
	if (nerr != 0) 
	{
		printf("Cannot access header\n");
		exit(nerr) ;
	}
	
	/* Begin iteration */
	for(n=0; n<STEP; n++)
	{
		h[n] = hmin+n*deltah;
		for(m=0; m<STEP; m++)
		{
			/* Calculation of Vs */
			k[m] = kmin+m*deltak;
			vs = vp/k[m];
			
			/* Predicted arrival times */
			vsterm = sqrtf(powf(vs,-2.0)-powf(p,2.0));
			vpterm = sqrtf(powf(vp,-2.0)-powf(p,2.0));
			tps = h[n]*(vsterm-vpterm);
			tppps = h[n]*(vsterm+vpterm);
			tppss = 2*h[n]*vsterm;
			
			/* Arrival times to data points */
			int dtps = start+tps*freq;
			int dtppps = start+tppps*freq;
			int dtppss = start+tppss*freq;
			
			/* Stack */
			s[n][m] = w1*array[dtps]+w2*array[dtppps]-w3*array[dtppss];
		}
	}	
	
	/* Write results to file */
	fout = fopen(output, "w");
	for(n=0; n<STEP; n++)
	{
		for(m=0; m<STEP; m++)
		{
			fprintf(fout, "%2.2f,%1.3f,%1.6f\n", h[n], k[m], s[n][m]);
		}
	}
	fclose(fout);;
	return 0;
}
