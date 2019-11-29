/********************************************************************
 *
 * HOPPING MODEL
 *
 * Author: R.H.J. Gerritsen
 *
 * Created on: 23-11-2019
 *
 *******************************************************************/

#include <time.h>
#include "main.h"
#include "functions.h"

double avgMobility;

void runExperiment(void)
{
	/************************************************************************************
	 * Procedure: Runs a Monte Carlo simulation for electron hopping
	 ************************************************************************************/

	/* Initialize the simulation and time the process*/
	clock_t tic = clock();
	initialise();

	clock_t toc1 = clock();
	printf("--- Setup Time: %f seconds ---\n", (1.0 * toc1 - 1.0 * tic) / CLOCKS_PER_SEC);

	/* Running the actual simulation */
	for (int step = 0; step < nrOfSteps; step++){

		/* Give some feedback on the progress */
		if ((step + 1) % (nrOfSteps / 100) == 0) {
			printf("\r progress: %3.0f percent", 100.0 * (step + 1) / (nrOfSteps));
			fflush(stdout);
		}

		updateNextEventList();

		executeEvent();
	}
	printf("\n");

	/* Generate output and free memory*/
	outputManager_mobility();
	free_memory();

	/* Print total duration to console*/
	clock_t toc2 = clock();
	printf("--- Simulation Time: %f seconds ---\n", (1.0 * toc2 - 1.0 * toc1) / CLOCKS_PER_SEC);
	printf("*************************************************************\n");
}

int main(void)
{
	/* We read the model parameters first ... */
	readModelParameters("../input/modelParameters.txt");
	
	clock_t toc_main;
	clock_t tic_main = clock();

	/* ... such that we can overwrite them if we please to run multiple experiments. */
	double sigma[] = {1.0, 2.0, 3.0, 5.0,};
	double fieldStrength[] = {0.00001, 0.00007, 0.00012};
	int numberOfParticleList[] = {1, 10,25,50, 100, 500, 1000}; 
	
	for(int k = 0; k < 4; k++) {
		for (int j = 0; j < 3; j++)	{
			for (int i = 0; i < 7; i++)	{
				DOS_sigma[0] = sigma[k]*kBT;
				E_Field = fieldStrength[j];
				nrOfParticles = numberOfParticleList[i];
				nrOfParticlesOfType[0] = numberOfParticleList[i];
				runExperiment(); // Actual run of the experiment
				toc_main = clock();
				printf("--- Total Time: %f seconds ---\n", (1.0 * toc_main- 1.0 * tic_main) / CLOCKS_PER_SEC);
				printf("*************************************************************\n");
			}
		}
	}

	free_memory_parameter_reading();

	return (0);
}