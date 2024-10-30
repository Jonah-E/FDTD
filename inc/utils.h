#ifndef __UTILS_H__
#define __UTILS_H__
#ifdef __cplusplus
extern "C" {
#endif
#include "options.h"
#include <stdbool.h>

enum time_categories {
  TOTAL_TIME = 0,
  CUDA_DIFF_TIME,
  LAST_TIME, // Placholder for calculation below, do not use
};

#define TOTAL_NR_TIMES (LAST_TIME - TOTAL_TIME)

void reset_times(double* time_elapsed);

/* Populate a given vector memory location with values from a given range.*/
void generateRandVector(DataType* v, int v_len, DataType min, DataType max,
                        unsigned int seed);

/*Get the current CPU time in seconds as a double.*/
double getCpuSeconds(void);

/* Print the header for the time data to std.*/
void print_header();

/* Print time data to std.*/
void print_times(const struct options* opt, double* times, unsigned int len);

DataType euclicianNormVector(DataType* vectorA, unsigned int length);
#ifdef __cplusplus
}
#endif
#endif /*__UTILS_H__*/
