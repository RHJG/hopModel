/********************************************************************
 *
 * HOPPING MODEL
 *
 * Author: R.H.J. Gerritsen
 *
 * Created on: 23-11-2019
 *
 *******************************************************************/

/* RNG SEED*/
extern long int SEED;

/* General Model Parameters */
extern int nrOfSteps;
extern int maxNeighbours;
extern int nrOfSites;
extern double rCutOff;
extern int nrOfParticles;
extern double Xmax;

/* Miller-Abrahams Rate */
extern double v0;
extern double alpha;
extern double kBT;

/* Density of States */
extern double DOS_mu;
extern double DOS_sigma;

/* Electric Field */
extern double E_Field;