#include "utils.h"

/* CUDA Kernel to multipy a vector with 1.25 for a number of iterations.*/
__global__ void vectorIterMult(DataType* v, unsigned int v_len,
                               unsigned int iter)
{
  const int idx = threadIdx.x + blockDim.x * blockIdx.x;

  if (idx >= v_len)
    return;

  for (unsigned int i = 0; i < iter; ++i) {
    v[idx] = 1.00005 * v[idx];
  }
}

