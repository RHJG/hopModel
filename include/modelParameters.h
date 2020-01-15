/********************************************************************
 *
 * HOPPING MODEL
 *
 * Author: R.H.J. Gerritsen
 *
 * Created on: 08-01-2020
 * 
 *******************************************************************/

/* RNG SEED*/
extern long int SEED;

/* General Model Parameters */
extern int nrOfSteps;
extern int maxNeighbours;
extern int nrOfSites;
extern double rCutOff;
extern int nrOfpTypes;
extern double* DOS_mu;
extern double* DOS_sigma;
extern int* nrOfParticlesOfType;
extern double* chargeOfType;
extern double* v0OfType;
extern double* alphaOfType;

extern int nrOfParticles; // computed by program does not need to be defined in modelParameters.txt
extern double maxXYZ[3];

/* Miller-Abrahams Rate */
extern double v0;
extern double alpha;
extern double kBT;

/* Electric Field */
extern double E_Field;
extern double electrodeEnergy_neg;