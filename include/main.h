/********************************************************************
 *
 * HOPPING MODEL
 *
 * Author: R.H.J. Gerritsen
 *
 * Created on: 08-01-2020
 * 
 *******************************************************************/

/* Generally used libraries */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "random.h"

/* User defined header files */
#include "modelParameters.h"

/* Particle Types */
#define ELEC 0
#define HOLE 1
#define SING 2
#define TRIP 3
#define CT   4
