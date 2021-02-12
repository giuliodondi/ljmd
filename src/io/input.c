#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "io.h"

#if defined(MPI_ENABLED)
#include "mpi_headers/mpi_comm.h"
#include "mpi_headers/mpi_utils.h"
#include <mpi.h>
#endif

/* helper function: read a line and then return
   the first string with whitespace stripped off */
int get_a_line(FILE *fp, char *buf) {
        char tmp[BLEN], *ptr;

        /* read a line and cut of comments and blanks */
        if (fgets(tmp, BLEN, fp)) {
                int i;

                ptr = strchr(tmp, '#');
                if (ptr)
                        *ptr = '\0';
                i = strlen(tmp);
                --i;
                while (isspace(tmp[i])) {
                        tmp[i] = '\0';
                        --i;
                }
                ptr = tmp;
                while (isspace(*ptr)) {
                        ++ptr;
                }
                i = strlen(ptr);
                strcpy(buf, tmp);
                return 0;
        } else {
                perror("problem reading input");
                return -1;
        }
        return 0;
}

/* 	handles parsing of initialisation file read through stdin
        allocates the position and velocity arrays and reads
        the restart file with the ICs
*/
int initialise(mdsys_t *sys, FILE *infile, file_names *fnames, int *nprint) {

        char line[BLEN], restfile[BLEN];
        FILE *fp;
        int i;

#if defined(MPI_ENABLED)
        double send_double_arr[8];
        if (sys->proc_id == 0) {
#endif

                /* read input file */
                if (get_a_line(infile, line))
                        return 1;
                sys->natoms = atoi(line);
                if (get_a_line(infile, line))
                        return 1;
                sys->mass = atof(line);
                if (get_a_line(infile, line))
                        return 1;
                sys->epsilon = atof(line);
                if (get_a_line(infile, line))
                        return 1;
                sys->sigma = atof(line);
                if (get_a_line(infile, line))
                        return 1;
                sys->rcut = atof(line);
                if (get_a_line(infile, line))
                        return 1;
                sys->box = atof(line);
                if (get_a_line(infile, restfile))
                        return 1;
                if (get_a_line(infile, fnames->trajfile))
                        return 1;
                if (get_a_line(infile, fnames->ergfile))
                        return 1;
                if (get_a_line(infile, line))
                        return 1;
                sys->nsteps = atoi(line);
                if (get_a_line(infile, line))
                        return 1;
                sys->dt = atof(line);
                if (get_a_line(infile, line))
                        return 1;
                *nprint = atoi(line);

#if defined(MPI_ENABLED)

                send_double_arr[0] = sys->natoms;
                send_double_arr[1] = sys->mass;
                send_double_arr[2] = sys->epsilon;
                send_double_arr[3] = sys->sigma;
                send_double_arr[4] = sys->rcut;
                send_double_arr[5] = sys->box;
                send_double_arr[6] = sys->nsteps;
                send_double_arr[7] = sys->dt;
        }

        /* now broadcast the struct to all processes */
        MPI_Bcast(send_double_arr, 8, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        /* copy all the data into the structs */
        if (sys->proc_id != 0) {
                sys->natoms = (int)send_double_arr[0];
                sys->mass = send_double_arr[1];
                sys->epsilon = send_double_arr[2];
                sys->sigma = send_double_arr[3];
                sys->rcut = send_double_arr[4];
                sys->box = send_double_arr[5];
                sys->nsteps = (int)send_double_arr[6];
                sys->dt = send_double_arr[7];
        }

        init_segments(sys->nprocs, sys->proc_id, sys->proc_seg, sys->natoms);
        printf("%d\n", sys->natoms);
        printf("%f\n", sys->rcut);
        printf("%d\n", sys->proc_seg->size);

#endif

        /* allocate memory */
        sys->rx = (double *)malloc(sys->natoms * sizeof(double));
        sys->ry = (double *)malloc(sys->natoms * sizeof(double));
        sys->rz = (double *)malloc(sys->natoms * sizeof(double));
		sys->fx = (double *)malloc(sys->natoms * sizeof(double));
        sys->fy = (double *)malloc(sys->natoms * sizeof(double));
        sys->fz = (double *)malloc(sys->natoms * sizeof(double));
#if !defined(MPI_ENABLED)
        sys->vx = (double *)malloc(sys->natoms * sizeof(double));
        sys->vy = (double *)malloc(sys->natoms * sizeof(double));
        sys->vz = (double *)malloc(sys->natoms * sizeof(double));
        
#else
        int seg_size = sys->proc_seg->size;
        if (sys->proc_id > 0) {
                sys->vx = (double *)malloc(seg_size * sizeof(double));
                sys->vy = (double *)malloc(seg_size * sizeof(double));
                sys->vz = (double *)malloc(seg_size * sizeof(double));

        } else {
                sys->vx = (double *)malloc(sys->natoms * sizeof(double));
                sys->vy = (double *)malloc(sys->natoms * sizeof(double));
                sys->vz = (double *)malloc(sys->natoms * sizeof(double));
#endif

        /* read restart */
        fp = fopen(restfile, "r");
        if (fp) {
                for (i = 0; i < sys->natoms; ++i) {
                        fscanf(fp, "%lf%lf%lf", sys->rx + i, sys->ry + i,
                               sys->rz + i);
                }
                for (i = 0; i < sys->natoms; ++i) {
                        fscanf(fp, "%lf%lf%lf", sys->vx + i, sys->vy + i,
                               sys->vz + i);
                }
                fclose(fp);
                // azzero(sys->fx, sys->natoms);
                // azzero(sys->fy, sys->natoms);
                // azzero(sys->fz, sys->natoms);
        } else {
                perror("cannot read restart file");
                return 3;
        }

#if defined(MPI_ENABLED)
}

// MPI_Barrier(MPI_COMM_WORLD);
mpi_send_pos_vel(sys);

#endif

return 0;
}