#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mpi_headers/mpi_utils.h"

void mpi_hello(int proc_id) { printf("hello from %d\n", proc_id); }

void split_dimension(const int nprocs, const int size, int splitting[3]) {

        // determine the two possible sizes of cell along the dimension
        float baselines = ((float)size) / ((float)nprocs);

        int h1 = (int)ceil(baselines);
        int h2 = (int)floor(baselines);

        int n1 = 0, n2 = 0, tot_size;

		int i;
        for (i = 0; i < nprocs; ++i) {
                n1 = (nprocs - i);
                n2 = i;
                tot_size = n1 * h1 + n2 * h2;
                if (tot_size == size) {
                        break;
                }
        }

        splitting[0] = n1;
        splitting[1] = h1;
        splitting[2] = h2;
}

void init_segments(const int nprocs, const int proc_id, arr_seg_t *proc_seg,
                   const int size) {

        int splitting[3];

        split_dimension(nprocs, size, splitting);

        if (proc_id < splitting[0]) {
                proc_seg->idx = proc_id * splitting[1];
                proc_seg->size = splitting[1];
        } else {
                proc_seg->idx = splitting[0] * splitting[1] +
                                (proc_id - splitting[0]) * splitting[2];
                proc_seg->size = splitting[2];
        }
        memcpy(proc_seg->splitting, splitting, 3 * sizeof(int));
}

extern void mpi_collective_comm_arrays(const int nprocs,
                                       const int *const splitting, int *count,
                                       int *offsets) {
		int p;
        for (p = 0; p < nprocs; ++p) {
                if (p < splitting[0]) {
                        count[p] = splitting[1];
                        offsets[p] = p * splitting[1];
                } else {
                        count[p] = splitting[2];
                        offsets[p] = splitting[0] * splitting[1] +
                                     (p - splitting[0]) * splitting[2];
                }
        }
}
