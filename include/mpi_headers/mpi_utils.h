#ifndef MPI_UTIL_H
#define MPI_UTIL_H

#ifdef __cplusplus
extern "C" {
#endif

struct arr_seg {
        int idx;
        int size;
        int *splitting;
};
typedef struct arr_seg arr_seg_t;

extern int right_triangle_area(int l);

extern void max_min_arr(int *arr, int length, int *max_out, int *min_out);

extern int max_min_index(int *arr, int length, int swch);

extern void split_triangle_equal_areas(int size, int nprocs, int *segment,
                                       int *segment_areas);

extern void init_segments(const int nprocs, const int proc_id,
                          arr_seg_t *proc_seg, const int size);

extern void mpi_collective_comm_arrays(const int nprocs,
                                       const int *const splitting, int *count,
                                       int *offset);

#ifdef __cplusplus
}
#endif

#endif