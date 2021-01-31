/*
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * baseline c version.
 */

#include <stdio.h>

// custom header files
#include "common.h"
#include "io.h"
#include "physics.h"

/* main */
int main(int argc, char **argv) {
        int nprint;
        file_names fnames;
        FILE *traj, *erg;
        mdsys_t sys;
        double t_start;

        printf("LJMD version %3.1f\n", LJMD_VERSION);

        t_start = wallclock();

        char init_err = initialise(&sys, stdin, &fnames, &nprint);
        if (init_err) {
                return init_err;
        }

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
                velverlet(&sys);
				ekin(&sys);
        }
        /**************************************************/

        /* clean up: close files, free memory */
        printf("Simulation Done. Run time: %10.3fs\n", wallclock() - t_start);
        fclose(erg);
        fclose(traj);

        cleanup_mdsys(&sys);

        return 0;
}
