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
#include "variables.h"
#include "functions.h"

int* int_1D_array(int np)
/* create an 1D array with size [np] of type int */
{
	int* a;

	a = (int*)calloc(np, sizeof(int));

	return a;
}

void free_int_1D_array(int* a) {
	/* free an 1D array with size [np] of type int */
	free(a);
}

double* double_1D_array(int np)
/* create an 1D array with size [np] of type double */
{
	double* a;

	a = (double*)calloc(np, sizeof(double));

	return a;
}

void free_double_1D_array(double* a) {
	/* free an 1D array with size [np] of type double */
	free(a);
}

int** int_2D_matrix(int nm, int np)
/* create an 2D matrix with size [nm, np] of type int */
{
	int i;
	int** m;

	m = (int**)calloc(nm, sizeof(int*));
	if (m) {
		for (i = 0; i < nm; i++)
			m[i] = (int*)calloc(np, sizeof(int));
	}
	return m;
}

void free_int_2D_matrix(int** m, int nm, int np)
/* frees an 2D matrix with size [nm, np] of type int */
{
	int i;

	for (i = 0; i < nm; i++)
		free(m[i]);

	free(m);
}

double** double_2D_matrix(int nm, int np)
/* create an 2D matrix with size [nm, np] of type double */
{
	int i;
	double** m;

	m = (double**)calloc(nm, sizeof(double*));
	if (m) {
		for (i = 0; i < nm; i++)
			m[i] = (double*)calloc(np, sizeof(double));
	}

	return m;
}

void free_double_2D_matrix(double** m, int nm, int np)
/* frees an 2D matrix with size [nm, np] of type int */
{
	int i;

	for (i = 0; i < nm; i++)
		free(m[i]);

	free(m);
}


void alloc_memory(void) {
	/* allocate memory for key-veriables*/

	/* Graph, Rates and Site Locations */
	neighbourList = int_2D_matrix(nrOfSites, maxNeighbours);
	rateList = double_2D_matrix(nrOfSites, maxNeighbours);
	siteXYZ = double_2D_matrix(nrOfSites, 3);
	siteEnergy = double_1D_array(nrOfSites);

	/* Next Event List */
	nextEventRates = double_1D_array(nrOfParticles * maxNeighbours);
	nextEvent = int_2D_matrix(nrOfParticles * maxNeighbours, 2);

	/* Store results */
	siteOccupation = double_1D_array(nrOfSites);
	lengthOfPath = double_2D_matrix(nrOfParticles, 3);
	driftVelocity = double_1D_array(nrOfParticles);
	mobility = double_1D_array(nrOfParticles);

	/* State variables */
	particlePositions = int_1D_array(nrOfParticles);
}


void free_memory(void) {
	/* Free the memory of the key-variables
	 * Note the this list should be equal to the memory allocation list.
	 *************************************************************************/
	free_int_2D_matrix(neighbourList, nrOfSites, maxNeighbours);
	free_double_2D_matrix(rateList, nrOfSites, maxNeighbours);
	free_double_2D_matrix(siteXYZ,nrOfSites,3);
	free_double_1D_array(nextEventRates);
	free_int_2D_matrix(nextEvent, nrOfParticles * maxNeighbours, 2);
	free_double_1D_array(siteEnergy);
	free_double_1D_array(siteOccupation);
	free_double_2D_matrix(lengthOfPath, nrOfParticles, 3);
	free_int_1D_array(particlePositions);
	free_double_1D_array(driftVelocity);
	free_double_1D_array(mobility);
}