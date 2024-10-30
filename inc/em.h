#ifndef __EM_H_
#define __EM_H_

#include "options.h"
#ifdef __cplusplus
extern "C" {
#endif
#include "utils.h"
#include <math.h>

typedef struct fields {
  /*The dimensions of the cavity*/
  int Nx, Ny, Nz;

  /*Number of cells.*/
  DataType Lx, Ly, Lz;

  /*The dimensions of each cell*/
  DataType dx, dy, dz;

  /* The magic timestep*/
  DataType dt;

  /*Fields*/
  DataType *Ex, *Ey, *Ez;
  DataType *Hx, *Hy, *Hz;

} fields_t;

#define EPS0 (8.8541878 * pow(10, -12))
#define MU0 (4 * pow(10, -7) * M_PI)
#define C0 299792458

#ifdef __cplusplus
}
#endif
#endif /*__EM_H_*/
