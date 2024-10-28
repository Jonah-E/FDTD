#ifndef __OPTIONS_H__
#define __OPTIONS_H__

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

#define DataType double

struct options {
  int Nx,Ny,Nz;
  DataType Lx,Ly,Lz;
  int timesteps;
  int it_batch_size;
  int sampling;
  unsigned int seed;
  bool run_graph;
  bool run_cpu;
  bool print_options;
  bool print_header;
};

int parse_arguments(struct options*, int, char**);
void print_options(const struct options*);

#ifdef __cplusplus
}
#endif
#endif /*__OPTIONS_H__*/
