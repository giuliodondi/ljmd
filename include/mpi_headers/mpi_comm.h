
#include "physics.h"

#ifdef __cplusplus
extern "C" {
#endif
	
extern void mpi_reduce_forces(mdsys_t *sys, const int *count,
                              const int *offsets);

extern void mpi_send_pos_vel(mdsys_t *sys);

extern void mpi_exchange_positions(mdsys_t *sys, const int *count,
                                   const int *offsets);

extern void mpi_reduce_UKT(mdsys_t *sys);

#ifdef __cplusplus
}
#endif