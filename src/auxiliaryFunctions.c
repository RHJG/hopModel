#include "main.h"

double hopRate(double dist, double distX, double charge, double energy1, double energy2, double v0, double alpha) {
	/* Compute the Miller-Abrahams hop rate, for given seperation and energies */
	if (energy2 - energy1 + E_Field*charge*(distX) <= 0) {
		return(v0 * exp(-2 * alpha * dist));
	}
	else {
		return(v0 * exp(-2 * alpha * dist - (energy2 - energy1 + E_Field * charge * (distX)) / kBT));
	}
}

double norm2(double dx, double dy, double dz) {
	/* returns the 2 norm of a 3 vector (x,y,z). */
	return(sqrt(dx * dx + dy * dy + dz * dz));
}

double distance_PBC_3D(double a[], double b[]) {
	/* compute the distance between 2 3-vectors a and b corrected for PBC. */
	double dx, dy, dz;
	dx = b[0] - a[0] - floor((b[0] - a[0]) / maxXYZ[0] + 0.5) * maxXYZ[0]; 
	dy = b[1] - a[1] - floor((b[1] - a[1]) / maxXYZ[1] + 0.5) * maxXYZ[1];
	dz = b[2] - a[2] - floor((b[2] - a[2]) / maxXYZ[2] + 0.5) * maxXYZ[2];
	return(norm2(dx, dy, dz));
}

double distance_PBC_3D_noX(double a[], double b[]) {
	/* compute the distance between 2 3-vectors a and b corrected for PBC. */
	double dx, dy, dz;
	dx = b[0] - a[0]; //We break the periodic boundary condition along the x-direction uncomment for reenabling.
	dy = b[1] - a[1] - floor((b[1] - a[1]) / maxXYZ[1] + 0.5) * maxXYZ[1];
	dz = b[2] - a[2] - floor((b[2] - a[2]) / maxXYZ[2] + 0.5) * maxXYZ[2];
	return(norm2(dx, dy, dz));
}

double correctedPBC_1D(double dx, double maxL) {
	return(dx - floor(dx / maxL + 0.5) * maxL);
}