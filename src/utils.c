#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>

/* Support function to generate a random number in a given range.
 * Based on the solution presented in the following forum thread:
 * http://ubuntuforums.org/showthread.php?t=1717717&p=10618266#post10618266*/
static DataType randfrom(DataType min, DataType max) {
  DataType range = (max - min);
  DataType div = ((DataType) RAND_MAX) / range;
  return min + (rand() / div);
}

/* Populate a given vector memory location with values from a given range.*/
void generateRandVector(DataType *v, int v_len, DataType min, DataType max,
                        unsigned int seed) {
  if (0 == seed)
    srand(time(NULL));
  else
    srand(seed);

  for (int i = 0; i < v_len; ++i) {
    v[i] = randfrom(min, max);
  }
}

/*Get the current CPU time in seconds as a double.*/
double getCpuSeconds(void) {
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return ((double)tp.tv_sec + (double)tp.tv_usec * 1.e-6);
}


const char *time_categories_string[] = {
  "total_time",
  "cuda_diff",
};

void print_header() {
  printf("tag,build,mode,timesteps,it_batch_size,Nx,Ny,Nz,Lx,Ly,Lz");
  for(unsigned int i = 0; i < TOTAL_NR_TIMES; ++i){
    printf(",%s", time_categories_string[i]);
  }
  printf("\n");
}

void print_times(const struct options *opt, double *times, unsigned int len) {
  printf("gpu,%s, ", BUILD_VERSION);

  // Print the options
  printf("%s, %d, %d, %d, %d, %d, %e, %e, %e",
          opt->run_graph ? "graph" : "nograph",
          opt->timesteps,
          opt->it_batch_size,
          opt->Nx,
          opt->Ny,
          opt->Nz,
          opt->Lx,
          opt->Ly,
          opt->Lz
         );

  for (unsigned int i = 0; i < len; ++i) {
    printf(", %lf", times[i]);
  }
  printf("\n");
}

void reset_times(double *time_elapsed) {
  for (int i = 0; i < TOTAL_NR_TIMES; ++i) {
    time_elapsed[i] = 0.0;
  }
}

DataType euclicianNormVector(DataType *vectorA, unsigned int length) {
  DataType diffEu = 0;

  for (int i = 0; i < length; ++i) {
    diffEu += vectorA[i] * vectorA[i];
  }
  diffEu = sqrt(diffEu);

  return diffEu;
}


