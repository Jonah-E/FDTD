#ifndef __HOST_H_
#define __HOST_H_

#ifdef __cplusplus
extern "C" {
#endif
#include "utils.h"
#include "em.h"

fields_t *host_setup(int Nx, int Ny, int Nz,
                     DataType Lx, DataType Ly, DataType Lz,
                     unsigned int seed);

void host_teardown(void);

int cpu_kernel_run(fields_t *h_fields, int time_steps, int sample_freq);

#ifdef __cplusplus
}
#endif
#endif /*__HOST_H_*/
