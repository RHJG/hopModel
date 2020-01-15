/********************************************************************
 *
 * KMC MODEL
 *
 * Author: R.H.J. Gerritsen
 *
 * Created on: 08-01-2020
 *
 * "random" contains all functions that are necessary for the
 * random number generation and drawing from different distributions.
 * it contains the mersenne twister to draw uniform random variables
 * and uses this to supply uniform, normal and exponential random variables.
 *******************************************************************/

/* Initialize the random number generator with a seed s. 
 * NOTE: this function must always be called once before use of any other function in random. */
extern void init_genrand(unsigned long s); 

/* Generates a random number from a uniform distribution on [0,1). */
extern double genrand_Real(void);                

/* Generates a random number from a normal distribution with mean m and standard deviation s, using a Box-Muller transform. */ 
extern double genrand_Normal(double m, double s); 

/* Generates a random number from an exponential distribution with rate parameter r. */
extern double genrand_Exp(double r);              

