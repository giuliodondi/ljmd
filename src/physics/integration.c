#include "physics.h"

#if defined(MPI_ENABLED)
#include "mpi_headers/mpi_utils.h"
#endif

// split the verlet algorithm into two functions

/* first part: propagate velocities by half and positions by full step  */
void verlet_1(mdsys_t *sys) {
        int r_idx;
		int i;
#if defined(MPI_ENABLED)
        r_idx = sys->proc_seg->idx;
        for (i = 0; i < sys->proc_seg->size; ++i) {
#else
        r_idx = 0;
        for (i = 0; i < sys->natoms; ++i) {
#endif
                sys->vx[i] += 0.5 * sys->dt / mvsq2e * sys->fx[i] / sys->mass;
                sys->vy[i] += 0.5 * sys->dt / mvsq2e * sys->fy[i] / sys->mass;
                sys->vz[i] += 0.5 * sys->dt / mvsq2e * sys->fz[i] / sys->mass;
                sys->rx[r_idx] += sys->dt * sys->vx[i];
                sys->ry[r_idx] += sys->dt * sys->vy[i];
                sys->rz[r_idx] += sys->dt * sys->vz[i];
                ++r_idx;
        }
}

/* second part: propagate velocities by another half step */
void verlet_2(mdsys_t *sys) {
		int i;
#if defined(MPI_ENABLED)
        for (i = 0; i < sys->proc_seg->size; ++i) {
#else
        for (i = 0; i < sys->natoms; ++i) {
#endif
                sys->vx[i] += 0.5 * sys->dt / mvsq2e * sys->fx[i] / sys->mass;
                sys->vy[i] += 0.5 * sys->dt / mvsq2e * sys->fy[i] / sys->mass;
                sys->vz[i] += 0.5 * sys->dt / mvsq2e * sys->fz[i] / sys->mass;
        }
}
