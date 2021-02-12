#include <math.h>

#include "physics.h"

#if defined(MPI_ENABLED)
#include "mpi_headers/mpi_utils.h"
#endif

/* helper function: apply minimum image convention */
static double pbc(double x, const double boxby2) {
        while (x > boxby2)
                x -= 2.0 * boxby2;
        while (x < -boxby2)
                x += 2.0 * boxby2;
        return x;
}

/* compute forces */
void force(mdsys_t *sys) {
        double r, ffac;
        double rx, ry, rz;
        int i, j;

        /* zero energy and forces */
        sys->epot = 0.0;
        double pot_energy;
        azzero(sys->fx, sys->natoms);
        azzero(sys->fy, sys->natoms);
        azzero(sys->fz, sys->natoms);

        double force_x, force_y, force_z;

#if defined(MPI_ENABLED)
		printf("from %d , %d , %d \n", sys->proc_id , sys->proc_seg->idx, sys->proc_seg->idx + sys->proc_seg->size );
        for (i = sys->proc_seg->idx;
             i < (sys->proc_seg->idx + sys->proc_seg->size); ++i) {
                for (j = i + 1; j < (sys->natoms); ++j) {
#else
        for (i = 0; i < (sys->natoms); ++i) {
                for (j = 0; j < (sys->natoms); ++j) {

                        /* particles have no interactions with themselves */
                        if (i == j)
                                continue;

#endif

                        /* get distance between particle i and j */
                        rx = pbc(sys->rx[i] - sys->rx[j], 0.5 * sys->box);
                        ry = pbc(sys->ry[i] - sys->ry[j], 0.5 * sys->box);
                        rz = pbc(sys->rz[i] - sys->rz[j], 0.5 * sys->box);
                        r = sqrt(rx * rx + ry * ry + rz * rz);

                        /* compute force and energy if within cutoff */
                        if (r < sys->rcut) {
                                ffac = -4.0 * sys->epsilon *
                                       (-12.0 * pow(sys->sigma / r, 12.0) / r +
                                        6 * pow(sys->sigma / r, 6.0) / r);

                                pot_energy = 4.0 * sys->epsilon *
                                             (pow(sys->sigma / r, 12.0) -
                                              pow(sys->sigma / r, 6.0));

#if !defined(MPI_ENABLED)
                                pot_energy /= 2;
#endif
                                sys->epot += pot_energy;

                                force_x = rx / r * ffac;
                                force_y = ry / r * ffac;
                                force_z = rz / r * ffac;

                                sys->fx[i] += force_x;
                                sys->fy[i] += force_y;
                                sys->fz[i] += force_z;
#if defined(MPI_ENABLED)
                                sys->fx[j] -= force_x;
                                sys->fy[j] -= force_y;
                                sys->fz[j] -= force_z;
#endif
                        }
                }
        }
}

