#include <stdlib.h>

#include "io.h"
#include "mpi_headers/mpi_utils.h"
#include <mpi.h>

void mpi_reduce_forces(mdsys_t *sys, const int *count, const int *offsets) {

        /* reduction in place, two different signatures for master and children
         */
        if (sys->proc_id == 0) {
                MPI_Reduce(MPI_IN_PLACE, sys->fx, sys->natoms, MPI_DOUBLE,
                           MPI_SUM, 0, MPI_COMM_WORLD);
                MPI_Reduce(MPI_IN_PLACE, sys->fy, sys->natoms, MPI_DOUBLE,
                           MPI_SUM, 0, MPI_COMM_WORLD);
                MPI_Reduce(MPI_IN_PLACE, sys->fz, sys->natoms, MPI_DOUBLE,
                           MPI_SUM, 0, MPI_COMM_WORLD);
        } else {
                MPI_Reduce(sys->fx, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM,
                           0, MPI_COMM_WORLD);
                MPI_Reduce(sys->fy, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM,
                           0, MPI_COMM_WORLD);
                MPI_Reduce(sys->fz, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM,
                           0, MPI_COMM_WORLD);
        }
}

void mpi_send_pos_vel(mdsys_t *sys) {

        /*	v{x,y,z}buf will act as the send buffer for velocities
           distribution for proc 0 they will be swapped with the current
           sys->v{x,y,z} array pointers and new smaller velocityes arrays
           allocated
        */

        double *vxbuf = NULL;
        double *vybuf = NULL;
        double *vzbuf = NULL;

        if (sys->proc_id == 0) {
                vxbuf = sys->vx;
                vybuf = sys->vy;
                vzbuf = sys->vz;

                sys->vx =
                    (double *)malloc(sys->proc_seg->size * sizeof(double));
                sys->vy =
                    (double *)malloc(sys->proc_seg->size * sizeof(double));
                sys->vz =
                    (double *)malloc(sys->proc_seg->size * sizeof(double));
        }

        /* 	initialise arrays for MPI Scatterv
                count contains the number of elements per segment scattered
                offset indicates the beginning of the segment in the master
           buffer
        */

        int count[sys->nprocs];
        int offsets[sys->nprocs];

        mpi_collective_comm_arrays(sys->nprocs, sys->proc_seg->splitting, count,
                                   offsets);

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Request rxreq, ryreq, rzreq, vxreq, vyreq, vzreq;

        /*scatter all the positions */
        MPI_Ibcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD, &rxreq);
        MPI_Ibcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD, &ryreq);
        MPI_Ibcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD, &rzreq);

        MPI_Wait(&rxreq, MPI_STATUS_IGNORE);
        MPI_Wait(&ryreq, MPI_STATUS_IGNORE);
        MPI_Wait(&rzreq, MPI_STATUS_IGNORE);

        /* scatter all the velocites*/

        MPI_Iscatterv(vxbuf, count, offsets, MPI_DOUBLE, sys->vx,
                      sys->proc_seg->size, MPI_DOUBLE, 0, MPI_COMM_WORLD,
                      &vxreq);
        MPI_Iscatterv(vybuf, count, offsets, MPI_DOUBLE, sys->vy,
                      sys->proc_seg->size, MPI_DOUBLE, 0, MPI_COMM_WORLD,
                      &vyreq);
        MPI_Iscatterv(vzbuf, count, offsets, MPI_DOUBLE, sys->vz,
                      sys->proc_seg->size, MPI_DOUBLE, 0, MPI_COMM_WORLD,
                      &vzreq);

        MPI_Wait(&vxreq, MPI_STATUS_IGNORE);
        MPI_Wait(&vyreq, MPI_STATUS_IGNORE);
        MPI_Wait(&vzreq, MPI_STATUS_IGNORE);

        if (sys->proc_id == 0) {
                free(vxbuf);
                free(vybuf);
                free(vzbuf);
        }
}

void mpi_exchange_positions(mdsys_t *sys, const int *count,
                            const int *offsets) {

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Request rxreq, ryreq, rzreq;

        /*	allgather all the positions
                        each process communicates to everyone its own segment
                        and all segments are ressembled into a full positions
           array
                */

        MPI_Iallgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, sys->rx, count,
                        offsets, MPI_DOUBLE, MPI_COMM_WORLD, &rxreq);

        MPI_Iallgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, sys->ry, count,
                        offsets, MPI_DOUBLE, MPI_COMM_WORLD, &ryreq);

        MPI_Iallgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, sys->rz, count,
                        offsets, MPI_DOUBLE, MPI_COMM_WORLD, &rzreq);

        MPI_Wait(&rxreq, MPI_STATUS_IGNORE);
        MPI_Wait(&ryreq, MPI_STATUS_IGNORE);
        MPI_Wait(&rzreq, MPI_STATUS_IGNORE);
}

void mpi_reduce_UKT(mdsys_t *sys) {

        MPI_Allreduce(MPI_IN_PLACE, &sys->epot, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &sys->ekin, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &sys->temp, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
}
