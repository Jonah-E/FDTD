#ifndef __UTILS_H__
#define __UTILS_H__
#ifdef __cplusplus
extern "C" {
#endif
#include <stdbool.h>
#include "options.h"


enum time_categories {
  TOTAL_TIME,
  CUDA_DIFF_TIME,
  GRAPH_CREATION,
  TOTAL_LAUNCH_COST,
  EXEC_TIME,
};

#define TOTAL_NR_TIMES (1 + EXEC_TIME - TOTAL_TIME)

void reset_times(double *time_elapsed);

/* Populate a given vector memory location with values from a given range.*/
void generateRandVector(DataType* v, int v_len, DataType min, DataType max, unsigned int seed);

/*Get the current CPU time in seconds as a double.*/
double getCpuSeconds(void);

/* Print time data to std.*/
void print_times(const struct options* opt, double* times, unsigned int len, DataType results);

DataType euclicianNormVector(DataType *vectorA, unsigned int length);
#ifdef __cplusplus
}
#endif
#endif /*__UTILS_H__*/
