/*
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * baseline c version.
 */

#include <stdio.h>
#include <unistd.h>

// custom header files
#include "common.h"
#include "io.h"
#include "physics.h"

#if defined(MPI_ENABLED)
#include "mpi_headers/mpi_io.h"
#include "mpi_headers/mpi_utils.h"
#include <mpi.h>
#endif

/* main */
int main(int argc, char **argv) {
        int nprint;
        file_names fnames;
        FILE *traj, *erg;
        mdsys_t sys;
        double t_start;

#if defined(MPI_ENABLED)
        int nprocs, proc_id;
        MPI_Init(&argc, &argv);

        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

        mpi_hello(proc_id);

        arr_seg_t proc_seg;

#endif

#if defined(MPI_ENABLED)
        if (proc_id == 0) {
#endif
                printf("LJMD version %3.1f\n", LJMD_VERSION);
#if defined(MPI_ENABLED)
        }
#endif

        t_start = wallclock();

#if defined(MPI_ENABLED)
        char init_err = mpi_initialise(nprocs, proc_id, &proc_seg, &sys, stdin,
                                       &fnames, &nprint);

#else
        char init_err = initialise(&sys, stdin, &fnames, &nprint);
#endif
        if (init_err) {
                return init_err;
        }

#if defined(MPI_ENABLED)
        /*
        for (int i = 0; i < nprocs; ++i) {
                if (proc_id == i) {
                        for (int j = proc_seg.idx;
                             j < (proc_seg.idx + proc_seg.size); ++j) {
                                printf("%d  %20.8f %20.8f %20.8f\n", j + 1,
                                       sys.rx[j], sys.ry[j], sys.rz[j]);
                        }
                }
                sleep(0.1);
                MPI_Barrier(MPI_COMM_WORLD);
        }
        for (int i = 0; i < nprocs; ++i) {
                if (proc_id == i) {
                        for (int j = 0; j < (proc_seg.size); ++j) {
                                printf("%d  %20.8f %20.8f %20.8f\n",
                                       j + proc_seg.idx + 1, sys.vx[j],
                                       sys.vy[j], sys.vz[j]);
                        }
                }
                sleep(0.1);
                MPI_Barrier(MPI_COMM_WORLD);
        }
                */

#endif

#if defined(MPI_ENABLED)
        cleanup_mdsys(&sys);
        MPI_Finalize();
        return 0;
#endif

        /* initialize forces and energies.*/
        sys.nfi = 0;
        force(&sys);
        ekin(&sys);

        erg = fopen(fnames.ergfile, "w");
        traj = fopen(fnames.trajfile, "w");

        printf("Startup time: %10.3fs\n", wallclock() - t_start);
        printf("Starting simulation with %d atoms for %d steps.\n", sys.natoms,
               sys.nsteps);
        printf("     NFI            TEMP            EKIN                 EPOT  "
               "      "
               "      ETOT\n");
        output(&sys, erg, traj);

        /* reset timer */
        t_start = wallclock();

        /**************************************************/
        /* main MD loop */
        for (sys.nfi = 1; sys.nfi <= sys.nsteps; ++sys.nfi) {

                /* write output, if requested */
                if ((sys.nfi % nprint) == 0)
                        output(&sys, erg, traj);

                /* propagate system and recompute energies */
                /* use the split versin of Verlet algorithm*/
                verlet_1(&sys);
                force(&sys);
                verlet_2(&sys);

                ekin(&sys);
        }
        /**************************************************/

        /* clean up: close files, free memory */
        printf("Simulation Done. Run time: %10.3fs\n", wallclock() - t_start);
        fclose(erg);
        fclose(traj);

        cleanup_mdsys(&sys);

#if defined(MPI_ENABLED)
        MPI_Finalize();
#endif

        return 0;
}
