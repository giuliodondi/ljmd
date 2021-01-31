#ifndef PHYSICS_H
#define PHYSICS_H

#include "common.h"

#ifdef __cplusplus
extern "C" {
#endif

/* a few physical constants */
static const double kboltz = 0.0019872067; /* boltzman constant in kcal/mol/K */
static const double mvsq2e = 2390.05736153349; /* m*v^2 in kcal/mol */

/* structure to hold the complete information
 * about the MD system */
struct _mdsys {
        int natoms, nfi, nsteps;
        double dt, mass, epsilon, sigma, box, rcut;
        double ekin, epot, temp;
        double *rx, *ry, *rz;
        double *vx, *vy, *vz;
        double *fx, *fy, *fz;
};
typedef struct _mdsys mdsys_t;

extern void cleanup_mdsys(mdsys_t *sys);
extern void ekin(mdsys_t *sys);
extern void force(mdsys_t *sys);
extern void velverlet(mdsys_t *sys);


#ifdef __cplusplus
}
#endif

#endif
