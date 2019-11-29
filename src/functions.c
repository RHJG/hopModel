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

double hopRate(double dist, double distX, double charge, double energy1, double energy2) {
	/* Compute the Miller-Abrahams hop rate, for given seperation and energies */
	if (energy2 - energy1 + E_Field*charge*(distX) <= 0) {
		return(v0 * exp(-2 * alpha * dist));
	}
	else {
		return(v0 * exp(-2 * alpha * dist - (energy2 - energy1 + E_Field * (distX)) / kBT));
	}
}

double norm2(double dx, double dy, double dz) {
	/* returns the 2 norm of a 3 vector (x,y,z). */
	return(sqrt(dx * dx + dy * dy + dz * dz));
}

double distance_PBC_3D(double a[], double b[]) {
	/* compute the distance between 2 3-vectors a and b corrected for PBC. */
	double dx, dy, dz;
	dx = b[0] - a[0] - floor((b[0] - a[0]) / maxXYZ[0] + 0.5) * maxXYZ[0]; 
	dy = b[1] - a[1] - floor((b[1] - a[1]) / maxXYZ[1] + 0.5) * maxXYZ[1];
	dz = b[2] - a[2] - floor((b[2] - a[2]) / maxXYZ[2] + 0.5) * maxXYZ[2];
	return(norm2(dx, dy, dz));
}

double distance_PBC_3D_noX(double a[], double b[]) {
	/* compute the distance between 2 3-vectors a and b corrected for PBC. */
	double dx, dy, dz;
	dx = b[0] - a[0]; //We break the periodic boundary condition along the x-direction uncomment for reenabling.
	dy = b[1] - a[1] - floor((b[1] - a[1]) / maxXYZ[1] + 0.5) * maxXYZ[1];
	dz = b[2] - a[2] - floor((b[2] - a[2]) / maxXYZ[2] + 0.5) * maxXYZ[2];
	return(norm2(dx, dy, dz));
}

double correctedPBC_1D(double dx, double maxL) {
	return(dx - floor(dx / maxL + 0.5) * maxL);
}

double random_Normal(double m, double s) {
	/************************************************************************************
	 * Input   : m the mean (mu) and s the standard deviation (sigma) of the normal distribution
	 * Output  : A random number ~Normal(mu,sigma)
	 * Procedure: Uses a Box-Muller transform method and a rng to compute a realisation of a normal
	 ************************************************************************************/
	double x1, x2, w, y1;
	static double y2;
	static int use_last = 0;

	if (use_last)		        /* use value from previous call */
	{
		y1 = y2;
		use_last = 0;
	}
	else
	{
		do {
			x1 = 2.0 * genrand_real2() - 1.0;
			x2 = 2.0 * genrand_real2() - 1.0;
			w = x1 * x1 + x2 * x2;
		} while (w >= 1.0);

		w = sqrt((-2.0 * log(w)) / w);
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}
	return(m + y1 * s);
}


/**************************************
 *        UPDATE EVENT LIST           *
***************************************/
void updateNextEventList(void) {
	int i, j;
	int counter = 0;
	int c_neighbour;

	/* scan through all current particles and their neighbouring sites for possible next events */
	for (i = 0; i < nrOfParticles; i++) {
		for (j = 0; j < maxNeighbours; j++) {
			c_neighbour = neighbourList[particlePositions[i]][j];
			if (c_neighbour == -1) {
				break;
			}
			if (siteIsOccupied[particleTypes[i]][c_neighbour] == 0) { // if it is occupied it should not be added to the next event list.
				nextEventRates[counter] = rateList[particleTypes[i]][particlePositions[i]][j];
				nextEvent[counter][0] = i; // particle i goes to ...
				nextEvent[counter][1] = c_neighbour; // ... position c_neighbour
				counter++;
			}
		}
	}
	/* fill the rest of the next event list with -1's */
	for (i = counter; i < nrOfParticles * maxNeighbours; i++) {
		nextEventRates[i] = 0;
		nextEvent[i][0] = -1;
		nextEvent[i][1] = -1;
	}

}


/**************************************
 *           EVENT EXECUTION          *
***************************************/
void executeEvent(void) {
	int i, j;
	int newPosition;
	int c_particle, c_position, c_type; // we use the c_ to mean current, the particle we are looking at at that instant.
	double rateSum = 0;
	double runSum;
	double jumpTime;
	double select;

	
	/* Compute total rates of all events */
	rateSum = 0;
	for (i = 0; i < nrOfParticles * maxNeighbours; i++) {
		rateSum += nextEventRates[i];	
	}

	/* Update occupation times */
	jumpTime = -(1.0 / rateSum) * log(genrand_real2());
	totalTime += jumpTime;
	for (i = 0; i < nrOfParticles; i++) {
		siteOccupation[particleTypes[i]][particlePositions[i]] += jumpTime;
	}

	/* select and execute the event */
	runSum = 0;
	select = rateSum * genrand_real2();
	for (i = 0; i < nrOfParticles * maxNeighbours; i++) {
		runSum += nextEventRates[i];
		if (runSum >= select) {
			if (nextEvent[i][1] == -1) {
				printf("-1 error in event selection.\n");
				exit(1);
			}
			/* Actual movement of particle */
			newPosition = nextEvent[i][1];
			c_particle = nextEvent[i][0];
			c_position = particlePositions[c_particle];
			c_type     = particleTypes[c_particle];
			siteIsOccupied[c_type][c_position] = 0; // free the old site
			siteIsOccupied[c_type][newPosition] = 1; // set the new position to occupied
			for (j = 0; j < 3; j++) { //update travelled distance 
				lengthOfPath[c_particle][j] += correctedPBC_1D(siteXYZ[newPosition][j] - siteXYZ[c_position][j], maxXYZ[j]);
			}
			particlePositions[c_particle] = newPosition;
			break;
		}
	}
}