/********************************************************************
 *
 * KMC MODEL
 *
 * Author: R.H.J. Gerritsen
 *
 * Created on: 08-01-2020
 * 
 * The outputManager handles all the output to files. 
 *******************************************************************/

#include "main.h"
#include "functions.h"
#include "variables.h"

void outputManager_mobility(void) {
	/* Generates output on the mobilities */
	char* filename = "../output/mobilities.txt";

	/* Open or create the output file */
	FILE* f = fopen(filename, "a");
	if (f == NULL)
	{
		printf("Error opening: %s!\n", filename);
		exit(1);
	}

	double driftSum = 0;
	double mobSum = 0;
	for (int i = 0; i < nrOfParticles; i++) {
		driftVelocity[i] = norm2(lengthOfPath[i][0] / totalTime, lengthOfPath[i][1] / totalTime, lengthOfPath[i][2] / totalTime);
		driftSum += driftVelocity[i];
		mobility[i] = driftVelocity[i] / E_Field;
		mobSum += mobility[i];
	}
	fprintf(f, "%4d %15E %15E %15E %15E %15E\n", nrOfParticles, alpha, E_Field, DOS_sigma[0]/kBT ,driftSum / nrOfParticles, mobSum / nrOfParticles);
	fclose(f);
}