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
#include <time.h>

#define max(a,b) (((a) > (b)) ? (a) : (b))

/*************************************
 *    MODEL PARAMETERS               *
 *************************************/
 /* RNG SEED*/
long int SEED;

/* General Model Parameters */
int nrOfSteps;
int maxNeighbours;
int nrOfSites;
double rCutOff;
int nrOfParticles;
double Xmax;

/* Miller-Abrahams Rate */
double v0;
double alpha;
double kBT;

/* Density of States */
double DOS_mu;
double DOS_sigma;

/* Electric Field */
double E_Field;

/*************************************
 *    GLOBAL VARIABLES               *
 *************************************/
/* Graph, Rates and Site Locations */
int** neighbourList;             // 2D-array containing the neighbourlist of the graph
double** rateList;               // 2D-array containing the rates corresponding to the neighbour in neighbourList
double** siteXYZ;                // 2D-array containing the (x,y,z) coordinate of every site
double* siteEnergy;              // 1D-array containing the siteEnergy

/* Cell list */
int* cellOfSite;                 // 1D-array containing the cell of a site
int** sitesInCell;               // 2D-array every row contains the sites in that cell

/* Next Event List */
double* nextEventRates;          // 1D-array row i corresponds to the rate of event i
int** nextEvent;                 // 2D-array row i gives the particle to update and to which site it is moved

/* Store results */
double* siteOccupation;           // 1D-array storing the site occupation
double** lengthOfPath; 
double* driftVelocity;
double* mobility;

/* State variables */
int* particlePositions;
double totalTime = 0;

int readSites(char* filename) {
	/************************************************************************************
	 * Input   : a filename
	 * Output  : Nothing important
	 * Procedure: Read the locations of the sites from a file
	 ************************************************************************************/
	FILE* fp;
	char str[1000];
	char c;
	float x, y, z;
	int count = 0;
	int nrElem;

	fp = fopen(filename, "r");
	if (fp == NULL) {
		printf("Could not open file %s", filename);
		exit(1);
	}

	/* Get number of sites (lines in file) and check against model parameters */
	for (c = getc(fp); c != EOF; c = getc(fp))
		if (c == '\n') // Increment count if this character is newline 
			count = count + 1;
	if (count != nrOfSites) {
		printf("Number of sites in file %s does not match nrOfSites in the model parameters.\n", filename);
		return 1;
	}
	else {
		printf("Number of sites matches.\n");
	}

	/* Now read the actual values */
	rewind(fp);
	for (count = 0; count < nrOfSites; count++) {
		fgets(str, 1000, fp);
		nrElem = sscanf(str, "%f %f %f", &x, &y, &z);
		siteXYZ[count][0] = x;
		siteXYZ[count][1] = y;
		siteXYZ[count][2] = z;
	}	

	fclose(fp);
	return 0;
}

int setupSiteEnergies(void) {
	/* Does what it says, assigns site energies to the siteEnergy variable. */
	for (int i = 0; i < nrOfSites; i++) {
		siteEnergy[i] = random_Normal(DOS_mu,DOS_sigma);
	}
	return(0);
}

int setupNeighboursAndRates(void) {
	/************************************************************************************
	 * Procedure: Create the neighbourlist and compute the jumping rates.
	 ************************************************************************************/
	double dist;
	int i, j, k;
	int c_neighbour;
	int maxN = 0;

	for (i = 0; i < nrOfSites; i++) {
		c_neighbour = 0;
		/* Add the neighbours whitin rCutOff to the neighbourlist. */
		for (j = 0; j < nrOfSites; j++) {
			if (i != j) {
				dist = distance_PBC_3D(siteXYZ[i],siteXYZ[j]);
				if (dist <= rCutOff) {
					neighbourList[i][c_neighbour] = j;
					rateList[i][c_neighbour] = hopRate(dist,  correctedPBC_1D(siteXYZ[j][0]-siteXYZ[i][0]), siteEnergy[i], siteEnergy[j]); 
					c_neighbour++;
					maxN = max(c_neighbour, maxN);
					if (c_neighbour == maxNeighbours) {
						printf("Number of maximum neighbours is too small!\n");
					}
				}
			}

		}
		/* fill the rest of the row with -1's */
		for (k = c_neighbour; k < maxNeighbours; k++) {
			neighbourList[i][k] = -1;
			rateList[i][k] = 0;
		}
	}
	printf("Maximum neighbours: %d\n", maxN);
	return(0);

}


void readModelParameters(char* filename) {
	/************************************************************************************
	* Input   : a filename
	* Procedure: Read the modelParamters from a file.
	************************************************************************************/
	FILE* fp;
	char str[1000];
	char junk[1000];
	char c;
	int count = 1;
	int nrElem;

	fp = fopen(filename, "r");
	if (fp == NULL) {
		printf("Could not open file %s", filename);
		exit(1);
	}

	/* Get number of sites (lines in file) and check against model parameters */
	for (c = getc(fp); c != EOF; c = getc(fp))
		if (c == '\n') // Increment count if this character is newline 
			count = count + 1;
	if (count != 13) {
		printf("Number of parameters in file %s does not match nr of parameters needed. found: %d\n", filename, count);
		exit(1);
	}

	/* Now read the actual values */
	rewind(fp);
	fgets(str, 1000, fp);
	nrElem = sscanf(str, "%s %d", &junk, &SEED);
	fgets(str, 1000, fp);
	nrElem = sscanf(str, "%s %d", &junk, &nrOfSteps);
	fgets(str, 1000, fp);
	nrElem = sscanf(str, "%s %d", &junk, &maxNeighbours);
	fgets(str, 1000, fp);
	nrElem = sscanf(str, "%s %d", &junk, &nrOfSites);
	fgets(str, 1000, fp);
	nrElem = sscanf(str, "%s %lf", &junk, &rCutOff);
	fgets(str, 1000, fp);
	nrElem = sscanf(str, "%s %d", &junk, &nrOfParticles);
	fgets(str, 1000, fp);
	nrElem = sscanf(str, "%s %lf", &junk, &Xmax);
	fgets(str, 1000, fp);
	nrElem = sscanf(str, "%s %lf", &junk, &v0);
	fgets(str, 1000, fp);
	nrElem = sscanf(str, "%s %lf", &junk, &alpha);
	fgets(str, 1000, fp);
	nrElem = sscanf(str, "%s %lf", &junk, &kBT);
	fgets(str, 1000, fp);
	nrElem = sscanf(str, "%s %lf", &junk, &DOS_mu);
	fgets(str, 1000, fp);
	nrElem = sscanf(str, "%s %lf", &junk, &DOS_sigma);
	fgets(str, 1000, fp);
	nrElem = sscanf(str, "%s %lf", &junk, &E_Field);
	fclose(fp);
}


void initialise(void) {
	/* Initialize all the variables. */
	int i, j;
	int select;

	totalTime = 0;
	
	alloc_memory(); 

	init_genrand(SEED); // Initialize random number generator
	
	readSites("../input/sites.txt"); // Read the site locations from a file

	setupSiteEnergies(); // Generate random site energies

	setupNeighboursAndRates(); // Setup the neigbourlist and compute all rates

	/* Initialize the result arrays with zeros. */
	for (i = 0; i < nrOfSites; i++) {
		siteOccupation[i] = 0;
	}
	for (i = 0; i < nrOfParticles; i++) {
		for (int j = 0; j < 3; j++) {
			lengthOfPath[i][j] = 0;
		}
	}
	/* Initialize the initial positions of the particles */
	printf("Initial particle positions: ");
	for (i = 0; i < nrOfParticles; i++) {
		select = (int) (1.0 * nrOfSites * genrand_real2());
		for (j = 0; j < i; j++) { // make sure all starting points are unique
			if (select == particlePositions[j]) {
				select = (int)(1.0 * nrOfSites * genrand_real2());
				j = 0;
			}
		}
		particlePositions[i] = select;
		printf("%d ", particlePositions[i]);
	}
	printf("\n");
	printf("Initialization and setup is done.\n");
}