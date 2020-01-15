/********************************************************************
 *
 * KMC MODEL
 *
 * Author: R.H.J. Gerritsen
 *
 * Created on: 08-01-2020
 *
 * The memoryManager handles all the memory allocation and freeing.
 *
 *******************************************************************/

#include "functions.h"
#include "main.h"
#include "variables.h"

/********************************************
 * GENERAL ARRAY AND MATRIX ALLOCATION
 * *****************************************/

int* int_1D_array(int np) {
    /* create an 1D array with size [np] of type int */
    int* a;

    a = (int*)calloc(np, sizeof(int));

    return a;
}

void free_int_1D_array(int* a) {
    /* free an 1D array with size [np] of type int */
    free(a);
}

double* double_1D_array(int np) {
    /* create an 1D array with size [np] of type double */
    double* a;

    a = (double*)calloc(np, sizeof(double));

    return a;
}

void free_double_1D_array(double* a) {
    /* free an 1D array with size [np] of type double */
    free(a);
}

int** int_2D_matrix(int nm, int np) {
    /* create an 2D matrix with size [nm, np] of type int */
    int i;
    int** m;

    m = (int**)calloc(nm, sizeof(int*));
    if (m) {
        for (i = 0; i < nm; i++) m[i] = (int*)calloc(np, sizeof(int));
    }
    return m;
}

void free_int_2D_matrix(int** m, int nm, int np) {
    /* frees an 2D matrix with size [nm, np] of type int */
    int i;

    for (i = 0; i < nm; i++) free(m[i]);

    free(m);
}

double** double_2D_matrix(int nm, int np) {
    /* create an 2D matrix with size [nm, np] of type double */
    int i;
    double** m;

    m = (double**)calloc(nm, sizeof(double*));
    if (m) {
        for (i = 0; i < nm; i++) m[i] = (double*)calloc(np, sizeof(double));
    }

    return m;
}

void free_double_2D_matrix(double** m, int nm, int np) {
    /* frees an 2D matrix with size [nm, np] of type int */
    int i;

    for (i = 0; i < nm; i++) free(m[i]);

    free(m);
}

int*** int_3D_matrix(int nx, int ny, int nz) {
    /* create a 3D matrix with size [nx, ny, nz] of type int */
    int i, j;
    int*** m;

    m = (int***)calloc(nx, sizeof(int**));
    if (m) {
        for (i = 0; i < nx; i++) {
            m[i] = (int**)calloc(ny, sizeof(int*));
        }
        for (i = 0; i < nx; i++) {
            for (j = 0; j < ny; j++) {
                m[i][j] = (int*)calloc(nz, sizeof(int));
            }
        }
    }
    return m;
}

void free_int_3D_matrix(int*** m, int nx, int ny, int nz) {
    /* free a 3D matrix with size [nx, ny, nz] of type int */
    int i, j;

    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            free(m[i][j]);
        }
    }

    for (i = 0; i < nx; i++) {
        free(m[i]);
    }

    free(m);
}

double*** double_3D_matrix(int nx, int ny, int nz) {
    /* create a 3D matrix with size [nx, ny, nz] of type double */
    int i, j;
    double*** m;

    m = (double***)calloc(nx, sizeof(double**));
    if (m) {
        for (i = 0; i < nx; i++) {
            m[i] = (double**)calloc(ny, sizeof(double*));
        }
        for (i = 0; i < nx; i++) {
            for (j = 0; j < ny; j++) {
                m[i][j] = (double*)calloc(nz, sizeof(double));
            }
        }
    }
    return (m);
}

void free_double_3D_matrix(double*** m, int nx, int ny, int nz) {
    /* free a 3D matrix with size [nx, ny, nz] of type int */
    int i, j;

    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            free(m[i][j]);
        }
    }

    for (i = 0; i < nx; i++) {
        free(m[i]);
    }

    free(m);
}

/*****************************************************
 * MEMORY ALLOCATION FOR THE VARIABLES AND PARAMETERS
 * ***************************************************/

void alloc_memory_parameter_reading(void) {
    DOS_mu = double_1D_array(nrOfpTypes);
    DOS_sigma = double_1D_array(nrOfpTypes);
    nrOfParticlesOfType = int_1D_array(nrOfpTypes);
    chargeOfType = double_1D_array(nrOfpTypes);
    v0OfType = double_1D_array(nrOfpTypes);
    alphaOfType = double_1D_array(nrOfpTypes);
}

void free_memory_parameter_reading(void) {
    free_double_1D_array(DOS_mu);
    free_double_1D_array(DOS_sigma);
    free_int_1D_array(nrOfParticlesOfType);
    free_double_1D_array(chargeOfType);
    free_double_1D_array(v0OfType);
    free_double_1D_array(alphaOfType);
}

void alloc_memory(void) {
    /* allocate memory for key-veriables*/

    /* Graph, Rates and Site Locations */
    neighbourList = int_2D_matrix(nrOfSites, maxNeighbours);
    rateList = double_3D_matrix(nrOfpTypes, nrOfSites, maxNeighbours);
    siteXYZ = double_2D_matrix(nrOfSites, 3);
    siteEnergy = double_2D_matrix(nrOfpTypes, nrOfSites);
    siteIsOccupied = int_2D_matrix(nrOfpTypes, nrOfSites);

    /* Next Event List */
    nextEventRates = double_1D_array(nrOfParticles * maxNeighbours);
    nextEvent = int_2D_matrix(nrOfParticles * maxNeighbours, 2);

    /* Store results */
    siteOccupation = double_2D_matrix(nrOfpTypes, nrOfSites);
    lengthOfPath = double_2D_matrix(nrOfParticles, 3);
    driftVelocity = double_1D_array(nrOfParticles);
    mobility = double_1D_array(nrOfParticles);

    /* State variables */
    particlePositions = int_1D_array(nrOfParticles);
    particleTypes = int_1D_array(nrOfParticles);
}

void free_memory(void) {
    /* Free the memory of the key-variables
     * Note that this list should be equal to the memory allocation list.
     *************************************************************************/

    /* Graph, Rates and Site Locations */
    free_int_2D_matrix(neighbourList, nrOfSites, maxNeighbours);
    free_double_3D_matrix(rateList, nrOfpTypes, nrOfSites, maxNeighbours);
    free_double_2D_matrix(siteXYZ, nrOfSites, 3);
    free_double_2D_matrix(siteEnergy, nrOfpTypes, nrOfSites);
    free_int_2D_matrix(siteIsOccupied, nrOfpTypes, nrOfSites);

    /* Next Event List */
    free_double_1D_array(nextEventRates);
    free_int_2D_matrix(nextEvent, nrOfParticles * maxNeighbours, 2);

    /* Store results */
    free_double_2D_matrix(siteOccupation, nrOfpTypes, nrOfSites);
    free_double_2D_matrix(lengthOfPath, nrOfParticles, 3);
    free_double_1D_array(driftVelocity);
    free_double_1D_array(mobility);

    /* State variables */
    free_int_1D_array(particlePositions);
    free_int_1D_array(particleTypes);
}