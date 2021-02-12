#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mpi_headers/mpi_utils.h"

/*calculates the area of a right triangle of cathetus l (int) */
int right_triangle_area(int l) { return l * (l + 1) / 2; }

/* finds the maximum and minimum value of an int array */
void max_min_arr(int *arr, int length, int *max_out, int *min_out) {

        int max = arr[0];
        int min = arr[0];
        for (int k = 1; k < length; ++k) {
                if (arr[k] < min) {
                        min = arr[k];
                } else if (arr[k] > max) {
                        max = arr[k];
                }
        }
        *max_out = max;
        *min_out = min;
}

/*
        finds the position of the maximum (or minimum) element in array
        swch=0 (and default value) for minimum, swch=1 for maximum
*/
int max_min_index(int *arr, int length, int swch) {
        int idx = 0;
        int prev_val = arr[0];
        if (swch == 1) {

                for (int i = 1; i < length; ++i) {
                        if (arr[i] > prev_val) {
                                prev_val = arr[i];
                                idx = i;
                        }
                }

        } else {

                for (int i = 1; i < length; ++i) {
                        if (arr[i] < prev_val) {
                                prev_val = arr[i];
                                idx = i;
                        }
                }
        }

        return idx;
}

/*
                        divide the matrix triangle in horizontal slices of equal
   area momentarily, the slices are calculated from bottom to top first
   calculate the "b segments" as the vertical length from the bottom right
   corner and the top of the k-th slice, they define sub-triangles
                        for instance an 8x8 upper triangle sliced in 3, sohwn is
   "b segment 1"
                        o o o o o o o o
                          o o o o o o o _
                            x x x x x x |
                                  x x x x x |
                                    # # # # |
                                          # # # |
                                            # # |
                                                  # |
                        calculate the triangle segments using three discrete
   approx methods calculate the areas of the resulting slice as well, not the
   area of the individual sub-triangles
                */

void split_triangle_equal_areas(int size, int nprocs, int *segment,
                                int *segment_areas) {

        int subtr_areas[3 * nprocs];

        for (int k = 0; k < nprocs; ++k) {
                double a_k = size * sqrt(((double)(k + 1)) / ((double)nprocs));
                segment[k] = floor(a_k);
                segment[nprocs + k] = ceil(a_k);
                segment[2 * nprocs + k] = round(a_k);
                subtr_areas[k] = right_triangle_area(segment[k]);
                subtr_areas[nprocs + k] =
                    right_triangle_area(segment[nprocs + k]);
                subtr_areas[2 * nprocs + k] =
                    right_triangle_area(segment[2 * nprocs + k]);

                segment_areas[k] = subtr_areas[k];
                segment_areas[nprocs + k] = subtr_areas[nprocs + k];
                segment_areas[2 * nprocs + k] = subtr_areas[2 * nprocs + k];
                if (k > 0) {
                        segment_areas[k] -= subtr_areas[k - 1];
                        segment_areas[nprocs + k] -=
                            subtr_areas[nprocs + k - 1];
                        segment_areas[2 * nprocs + k] -=
                            subtr_areas[2 * nprocs + k - 1];
                }
        }
}

void init_segments(const int nprocs, const int proc_id, arr_seg_t *proc_seg,
                   const int size) {

        int segment[3 * nprocs];
        int segment_areas[3 * nprocs];

        split_triangle_equal_areas(size, nprocs, segment, segment_areas);

        /*
                to determine the best approx method we calculate the delta area
           between the largest area slice and the smallest area slice since this
           determines process imbalancing choose the approx method that
           minimises this
        */

        int delta_area[3];
        for (int i = 0; i < 3; ++i) {
                int max_area, min_area;

                max_min_arr(&segment_areas[i * nprocs], nprocs, &max_area,
                            &min_area);

                delta_area[i] = max_area - min_area;
        }

        int best_method = max_min_index(delta_area, 3, 0);

        int *best_segment = &segment[nprocs * best_method];

        /*now calculate the actual slices heights and store them in the
         * splitting array */

        proc_seg->splitting[nprocs - 1] = best_segment[0];
        for (int k = 1; k < nprocs; ++k) {
                proc_seg->splitting[nprocs - k - 1] =
                    segment[k] - segment[k - 1];
        }

        proc_seg->size = proc_seg->splitting[proc_id];
        proc_seg->idx = size - best_segment[nprocs - proc_id - 1];
}

void mpi_collective_comm_arrays(const int nprocs,
                                       const int *const splitting, int *count,
                                       int *offsets) {

        memcpy(count, splitting, nprocs * sizeof(int));

        offsets[0] = 0;
        for (int p = 0; p < nprocs - 1; ++p) {
                offsets[p + 1] = offsets[p] + splitting[p];
        }
}