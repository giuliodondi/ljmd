#include "physics.h"

// split the verlet algorithm into two functions

/* first part: propagate velocities by half and positions by full step  */
void verlet_1(mdsys_t *sys) {
		int i;
        for (i = 0; i < sys->natoms; ++i) {
                sys->vx[i] += 0.5 * sys->dt / mvsq2e * sys->fx[i] / sys->mass;
                sys->vy[i] += 0.5 * sys->dt / mvsq2e * sys->fy[i] / sys->mass;
                sys->vz[i] += 0.5 * sys->dt / mvsq2e * sys->fz[i] / sys->mass;
                sys->rx[i] += sys->dt * sys->vx[i];
                sys->ry[i] += sys->dt * sys->vy[i];
                sys->rz[i] += sys->dt * sys->vz[i];
        }
}

/* second part: propagate velocities by another half step */
void verlet_2(mdsys_t *sys) {
		int i;
        for (i = 0; i < sys->natoms; ++i) {
                sys->vx[i] += 0.5 * sys->dt / mvsq2e * sys->fx[i] / sys->mass;
                sys->vy[i] += 0.5 * sys->dt / mvsq2e * sys->fy[i] / sys->mass;
                sys->vz[i] += 0.5 * sys->dt / mvsq2e * sys->fz[i] / sys->mass;
        }
}
