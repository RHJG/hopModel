/********************************************************************
 *
 * KMC MODEL
 *
 * Author: R.H.J. Gerritsen
 *
 * Created on: 08-01-2020
 *
 * "functions" contains all functions of the KMC model.
 *******************************************************************/

/* INITILIZATION */

extern void readModelParameters(char* filename);
extern void initialise(void);

/* EVENT HANDLING */

extern void updateNextEventList(void);
extern void executeEvent(void);

/* OUTPUT MANAGER */

/* Outputs mobilities to a file */
extern void outputManager_mobility(void); 

/* MEMORY MANAGER */

/* Allocate the memory necessary for the parameters */
extern void alloc_memory_parameter_reading(void);
/* Free the memory necessary for the parameters */
extern void free_memory_parameter_reading(void);
/* Allocate the memory for all key-variables (ratelists, nb lists etc.) */
extern void alloc_memory(void);
/* Free the memory for all key-variables (ratelists, nb lists etc.) */
extern void free_memory(void);


/* AUXILIARY FUNCTIONS */
extern double hopRate(double dist, double distX, double charge, double energy1, double energy2, double v0, double alpha);
extern double norm2(double dx, double dy, double dz);
extern double distance_PBC_3D(double a[], double b[]);
extern double distance_PBC_3D_noX(double a[], double b[]);
extern double correctedPBC_1D(double dx, double maxL);