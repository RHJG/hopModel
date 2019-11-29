/********************************************************************
 *
 * HOPPING MODEL
 *
 * Author: R.H.J. Gerritsen
 *
 * Created on: 23-11-2019
 *
 *******************************************************************/


/* Graph, Rates, Energies and Site Locations */
extern int** neighbourList;       // 2D-array containing the neighbourlist of the graph
extern double*** rateList;         // 3D-array containing the rates corresponding to the neighbour in neighbourList for a certain particle type: rateList[pType][currentPosition][indexOfNeighbourInNeighbourList]
extern double** siteXYZ;          // 2D-array containing the (x,y,z) coordinate of every site
extern double** siteEnergy;        // 2D-array containing the siteEnergy siteEnergy[pType][siteID=pos]
extern int** siteIsOccupied;       // 2D-array storing the current state of a site, occupied or not for different particle types siteIsOccupied[pType][pos=siteID]

/* Next Event List */
extern double* nextEventRates;       // 1D-array row i corresponds to the rate of event i
extern int** nextEvent;              // 2D-array row i gives the particle to update and to which site it is moved

/* Store results */
extern double** siteOccupation;    // 2D-array Contains the total time a site was occupied by particle type siteOccupation[pType][pos=siteID]
extern double** lengthOfPath;     // contains for every particle the total length of the path it travelled as a 3-vector 
extern double* driftVelocity;     // stores the driftvelocities of the particles
extern double* mobility;          // store the mobilities ..

/* State variables */
extern int* particlePositions;  // 1D-array containing the current sites of the particles particlePositions[pID]
extern int* particleTypes;      // 1D-array containing the particletype of particle i: particleTypes[pID]
extern double totalTime;        // is what is says

