/********************************************************************
 *
 * KMC MODEL
 *
 * Author: R.H.J. Gerritsen
 *
 * Created on: 08-01-2020
 *
 *
 *******************************************************************/

#include "main.h"
#include <time.h>
#include "functions.h"

void runExperiment(void) {

    /* Initialize the simulation and time the process*/
    clock_t tic = clock();
    initialise();
    clock_t toc1 = clock();
    printf("--- Setup Time: %f seconds ---\n",
           (1.0 * toc1 - 1.0 * tic) / CLOCKS_PER_SEC);

    /* Running the actual simulation */
    for (int step = 0; step < nrOfSteps; step++) {
        updateNextEventList();

        executeEvent();

        /* Give some feedback on the progress */
        if ((step + 1) % (nrOfSteps / 100) == 0) {
            printf("\r progress: %3.0f percent",
                   100.0 * (step + 1) / (nrOfSteps));
            fflush(stdout);
        }
    }
    printf("\n");

    /* Generate output and free memory*/
    outputManager_mobility();
    free_memory();

    /* Print total duration to console*/
    clock_t toc2 = clock();
    printf("--- Simulation Time: %f seconds ---\n",
           (1.0 * toc2 - 1.0 * toc1) / CLOCKS_PER_SEC);
    printf("*************************************************************\n");
}

int main(void) {

    /* We read the model parameters first ... */
    readModelParameters("../input/modelParameters.txt");
    /* ... such that we can overwrite them if we please to run multiple experiments. */

    /* Start timing */
    clock_t toc_main;
    clock_t tic_main = clock();

    /* Run the experiment */
    runExperiment();

    /* Report simulation duration */
    toc_main = clock();
    printf("--- Total Time: %f seconds ---\n",
           (1.0 * toc_main - 1.0 * tic_main) / CLOCKS_PER_SEC);
    printf("*************************************************************\n");

    /* Free last bits of memory */
    free_memory_parameter_reading();

    return (0);
}