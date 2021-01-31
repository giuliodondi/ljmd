#include "physics.h"

// split the verlet algorithm into two functions

/* first part: propagate velocities by half and positions by full step  */
static void verlet_1(mdsys_t *sys) {
        for (int i = 0; i < sys->natoms; ++i) {
                sys->vx[i] += 0.5 * sys->dt / mvsq2e * sys->fx[i] / sys->mass;
                sys->vy[i] += 0.5 * sys->dt / mvsq2e * sys->fy[i] / sys->mass;
                sys->vz[i] += 0.5 * sys->dt / mvsq2e * sys->fz[i] / sys->mass;
                sys->rx[i] += sys->dt * sys->vx[i];
                sys->ry[i] += sys->dt * sys->vy[i];
                sys->rz[i] += sys->dt * sys->vz[i];
        }
}

/* second part: propagate velocities by another half step */
static void verlet_2(mdsys_t *sys) {
        for (int i = 0; i < sys->natoms; ++i) {
                sys->vx[i] += 0.5 * sys->dt / mvsq2e * sys->fx[i] / sys->mass;
                sys->vy[i] += 0.5 * sys->dt / mvsq2e * sys->fy[i] / sys->mass;
                sys->vz[i] += 0.5 * sys->dt / mvsq2e * sys->fz[i] / sys->mass;
        }
}

/* velocity verlet */
void velverlet(mdsys_t *sys) {
        verlet_1(sys);
        force(sys);
        verlet_2(sys);
}