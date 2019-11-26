/********************************************************************
 *
 * HOPPING MODEL
 *
 * Author: R.H.J. Gerritsen
 *
 * Created on: 23-11-2019
 *
 *******************************************************************/

#include "main.h"
#include "functions.h"
#include "variables.h"

void outputManager_singleRun(void) {
	/* Energy versus occupation results */
	FILE* f = fopen("../output/energyVSocc.txt", "w");
	if (f == NULL)
	{
		printf("Error opening file energyVSocc file!\n");
		exit(1);
	}

	for (int i = 0; i < nrOfSites; i++) {
		fprintf(f, "%f %f\n", siteEnergy[i], siteOccupation[i]);
	}

	fclose(f);

	/* Electron mobility and drift velocities */
	f = fopen("mobilities.txt", "w");
	if (f == NULL)
	{
		printf("Error opening mobility file!\n");
		exit(1);
	}

	double driftSum=0;
	double mobSum=0;
	printf("%3s %15s %15s\n", "ID", "Drift vel.", "Mob.");
	for (int i = 0; i < nrOfParticles; i++) {
		driftVelocity[i] = norm2(lengthOfPath[i][0] / totalTime, lengthOfPath[i][1] / totalTime, lengthOfPath[i][2] / totalTime);
		driftSum += driftVelocity[i];
		mobility[i] = driftVelocity[i] / E_Field;
		mobSum += mobility[i];
		printf("%3d %15.4f %15.4f\n", i, driftVelocity[i], mobility[i]);
		fprintf(f, "%3d %15.6f %15.6f\n", i, driftVelocity[i], mobility[i]);
	}
	printf("%3s %15.4f %15.4f\n", "avg", driftSum / nrOfParticles, mobSum / nrOfParticles);
	fprintf(f, "%3s %15.6f %15.6f\n", "avg", driftSum / nrOfParticles, mobSum / nrOfParticles);

	fclose(f);
}

void outputManager_experiment(void) {
	FILE* f = fopen("../output/densityExperiment.txt", "a");
	if (f == NULL)
	{
		printf("Error opening file energyVSocc file!\n");
		exit(1);
	}

	double driftSum = 0;
	double mobSum = 0;
	printf("%3s %15s %15s\n", "ID", "Drift vel.", "Mob.");
	for (int i = 0; i < nrOfParticles; i++) {
		driftVelocity[i] = norm2(lengthOfPath[i][0] / totalTime, lengthOfPath[i][1] / totalTime, lengthOfPath[i][2] / totalTime);
		driftSum += driftVelocity[i];
		mobility[i] = driftVelocity[i] / E_Field;
		mobSum += mobility[i];
		printf("%3d %15.4f %15.4f\n", i, driftVelocity[i], mobility[i]);
	}
	printf("%3s %15.4f %15.4f\n", "avg", driftSum / nrOfParticles, mobSum / nrOfParticles);
	fprintf(f, "%15E %15E\n", E_Field, driftSum / nrOfParticles);
	avgMobility = driftSum / nrOfParticles;


	fclose(f);
}

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
	fprintf(f, "%4d %15E %15E %15E %15E\n", nrOfParticles, E_Field, DOS_sigma/kBT ,driftSum / nrOfParticles, mobSum / nrOfParticles);
	fclose(f);
}
