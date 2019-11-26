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
extern double** rateList;         // 2D-array containing the rates corresponding to the neighbour in neighbourList
extern double** siteXYZ;          // 2D-array containing the (x,y,z) coordinate of every site
extern double* siteEnergy;        // 1D-array containing the siteEnergy

/* Cell list */
extern int* cellOfSite;                // 1D-array containing the cell of a site
extern int** sitesInCell;              // 2D-array every row contains the sites in that cell

/* Next Event List */
extern double* nextEventRates;       // 1D-array row i corresponds to the rate of event i
extern int** nextEvent;              // 2D-array row i gives the particle to update and to which site it is moved

/* Store results */
extern double* siteOccupation;    // 1D_array Contains the total time a site was occupied
extern double** lengthOfPath;     // contains for every particle the total length of the path it travelled as a 3-vector 
extern double* driftVelocity;     // stores the driftvelocities of the particles
extern double* mobility;          // store the mobilities ..
extern double avgMobility;        // stores the average mobility

/* State variables */
extern int* particlePositions;  // 1D-array containing the current sites of the particles
extern double totalTime;        // is what is says

