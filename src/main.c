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

void runSingleSimulation(void) {
	/************************************************************************************
	 * Procedure: Runs a Monte Carlo simulation for electron hopping
	 ************************************************************************************/

	/* Initialize the simulation and time the process*/
	clock_t tic = clock();
	initialise();
	clock_t toc1 = clock();
	printf("--- Setup Time: %f seconds ---\n", (1.0 * toc1 - 1.0 * tic) / CLOCKS_PER_SEC);

	/* Running the actual simulation */
	for (int step = 0; step < nrOfSteps; step++) {

		/* Give some feedback on the progress */
		if ( (step+1) % (nrOfSteps / 100) == 0) printf("\r progress: %3.0f percent", 100.0 * (step+1) / (nrOfSteps));

		updateNextEventList();

		executeEvent();
	}

	/* Generate output and free memory*/
	outputManager_singleRun();
	free_memory();

	/* Print total duration to console*/
	clock_t toc2 = clock();
	printf("--- Simulation Time: %f seconds ---\n", (1.0 * toc2 - 1.0 * toc1) / CLOCKS_PER_SEC);
}

void runExperiment(void) {
	/************************************************************************************
	 * Procedure: Runs a Monte Carlo simulation for electron hopping
	 ************************************************************************************/

	 /* Initialize the simulation and time the process*/
	clock_t tic = clock();
	initialise();
	clock_t toc1 = clock();
	printf("--- Setup Time: %f seconds ---\n", (1.0 * toc1 - 1.0 * tic) / CLOCKS_PER_SEC);

	/* Running the actual simulation */
	for (int step = 0; step < nrOfSteps; step++) {

		/* Give some feedback on the progress */
		if ( (step+1) % (nrOfSteps / 100) == 0) {
			printf("\r progress: %3.0f percent", 100.0 * (step+1) / (nrOfSteps));
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


int main(void) {
	readModelParameters("../input/modelParameters.txt");

	int numberOfParticleList[] = { 1,5,10,25,50,100};

	for (int i = 0; i < 6; i++) {
		nrOfParticles = numberOfParticleList[i];
		runExperiment();
	}

	return(0);
}