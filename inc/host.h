#ifndef __HOST_H_
#define __HOST_H_

#ifdef __cplusplus
extern "C" {
#endif
#include "em.h"

fields_t* host_setup(int Nx, int Ny, int Nz, DataType Lx, DataType Ly,
                     DataType Lz, unsigned int seed, bool generate);

void host_teardown(fields_t* h_fields);

int cpu_kernel_run(fields_t* h_fields, int time_steps, int sample_freq);

bool compaire_fields(fields_t* field_1, fields_t* field_2, DataType tolerance);

#ifdef __cplusplus
}
#endif
#endif /*__HOST_H_*/
