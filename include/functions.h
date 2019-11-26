/********************************************************************
 *
 * HOPPING MODEL
 *
 * Author: R.H.J. Gerritsen
 *
 * Created on: 23-11-2019
 *
 *******************************************************************/

/* FROM THE INITILIZE.c FILE */
extern void readModelParameters(char* filename);
extern void initialise(void);

/* FROM FUNCTIONS.c */

 /* MAIN */
extern void updateNextEventList(void);
extern void executeEvent(void);
extern void outputManager_singleRun(void);
extern void outputManager_experiment(void);
extern void outputManager_mobility(void);

/* Memory Managing */
extern void alloc_memory(void);
extern void free_memory(void);

/* Auxiliary Functions */
extern double hopRate(double dist, double distX, double energy1, double energy2);
extern double norm2(double dx, double dy, double dz);
extern double distance_PBC_3D(double a[], double b[]);
extern double correctedPBC_1D(double dx);
extern double random_Normal(double m, double s);

/* Random Numbers (from mersenne twister) */
extern void init_genrand(unsigned long s); // set seed for mersenne twister
extern double genrand_real2(void); // generate random number on [0,1)
 
/* Diagnostics */
extern void printHeadRates(void);
extern void printHeadNeighbours(void);
extern void printAllNeighbours(int p);
extern void printNextEventRates(void);
extern void printAllRates(int p);

