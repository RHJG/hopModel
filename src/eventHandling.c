/********************************************************************
 *
 * KMC MODEL
 *
 * Author: R.H.J. Gerritsen
 *
 * Created on: 08-01-2020
 *
 * "eventHandling" contains all function to handle events: 
 * the setup (setting up the next event list) and execution.
 *******************************************************************/

#include "main.h"
#include "functions.h"
#include "variables.h"

/**************************************
 *        UPDATE EVENT LIST           
***************************************/
void updateNextEventList(void) {
	int i, j;
	int counter = 0;
	int c_neighbour;

	/* scan through all current particles and their neighbouring sites for possible next events */
	for (i = 0; i < nrOfParticles; i++) {
		if(particlePositions[i] >=0 ){
			for (j = 0; j < maxNeighbours; j++) {
				c_neighbour = neighbourList[particlePositions[i]][j];
				if (c_neighbour == -11) {
					break;
				}
				if (( c_neighbour >=0)){
					if (siteIsOccupied[particleTypes[i]][c_neighbour] == 0) { // if it is occupied it should not be added to the next event list.
						nextEventRates[counter] = rateList[particleTypes[i]][particlePositions[i]][j];
						nextEvent[counter][0] = i; // particle i goes to ...
						nextEvent[counter][1] = c_neighbour; // ... position c_neighbour
						counter++;
					}
				} else {
					nextEventRates[counter] = rateList[particleTypes[i]][particlePositions[i]][j];
					nextEvent[counter][0] = i; // particle i goes to ...
					nextEvent[counter][1] = c_neighbour; // ... position c_neighbour
					counter++;
				}
			}
		} 
	}
	/* fill the rest of the next event list with -11's */
	for (i = counter; i < nrOfParticles * maxNeighbours; i++) {
		nextEventRates[i] = 0;
		nextEvent[i][0] = -11;
		nextEvent[i][1] = -11;
	}

}

/**************************************
 *           EVENT EXECUTION          
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
	jumpTime = -(1.0 / rateSum) * log(genrand_Real());
	totalTime += jumpTime;
	for (i = 0; i < nrOfParticles; i++) {
		if(particlePositions[i] >=0){
			siteOccupation[particleTypes[i]][particlePositions[i]] += jumpTime;
		}
	}

	/* select and execute the event */
	runSum = 0;
	select = rateSum * genrand_Real();
	for (i = 0; i < nrOfParticles * maxNeighbours; i++) {
		runSum += nextEventRates[i];
		if (runSum >= select) {
			if (nextEvent[i][1] == -11) {
				printf(" \nNo next events, all particles at electrode!\n");
				exit(0);
			}
			/* Actual movement of particle */
			newPosition = nextEvent[i][1];
			c_particle = nextEvent[i][0];
			c_position = particlePositions[c_particle];
			c_type     = particleTypes[c_particle];
			siteIsOccupied[c_type][c_position] = 0; // free the old site
			if (newPosition>=0 ){
				siteIsOccupied[c_type][newPosition] = 1; // set the new position to occupied
				for (j = 0; j < 3; j++) { //update travelled distance 
					lengthOfPath[c_particle][j] += correctedPBC_1D(siteXYZ[newPosition][j] - siteXYZ[c_position][j], maxXYZ[j]);
				}
			}
			particlePositions[c_particle] = newPosition;
			break;
		}
	}
}