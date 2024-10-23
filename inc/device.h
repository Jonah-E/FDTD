#ifndef __DEVICE_H__
#define __DEVICE_H__

#ifdef __cplusplus
extern "C" {
#endif
#include "options.h"
#include "em.h"

fields_t *device_setup(fields_t *h_fields);

int device_kernel_run(fields_t *d_fields, int timesteps);

void device_get_fields(fields_t *h_fields, fields_t *d_fields);

void device_teardown(void);

#ifdef __cplusplus
}
#endif
#endif /*__DEVICE_H__*/
