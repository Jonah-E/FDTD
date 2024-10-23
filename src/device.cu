#include "device.h"
#include "host.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>

#if defined(__HIP)
#include "hip/hip_runtime.h"
#define cudax hip##x
#else
#define cudax cuda##x
#endif

#define printCudaError(cuda_returned_error_code)                               \
  { accErrorPrint((cuda_returned_error_code), __FILE__, __LINE__); }

inline void accErrorPrint(cudaError_t code, const char *file, int line) {
  fprintf(stderr, "ACC Error: %s (%d) %s %d\n", cudaGetErrorString(code), code,
          file, line);
}


#include "simple-kernel.cu"

static DataType *gd_vector;
static DataType *device_setup(DataType *h_vector, unsigned int v_len) {
  cudaError_t device_error;
  device_error = cudaMalloc(&gd_vector, sizeof(DataType) * v_len);
  if (device_error != cudaSuccess) {
    printCudaError(device_error);
    return NULL;
  }

  device_error = cudaMemcpy(gd_vector, h_vector, sizeof(DataType) * v_len,
                            cudaMemcpyHostToDevice);
  if (device_error != cudaSuccess) {
    printCudaError(device_error);
    return NULL;
  }

  return gd_vector;
}

void device_teardown(void) { cudaFree(gd_vector); }


#define TPB 1024
int device_kernel_run(const struct options *opt, DataType *d_vector) {
  cudaError_t device_error;

  dim3 block(TPB);
  dim3 grid((opt->number_of_threads + TPB - 1) / TPB);

#ifdef TIME_DETAILED
  double time_start;
  time_start = getCpuSeconds();
#endif
  for (unsigned int i = 0; i < opt->outer_iterations; ++i) {
    for (unsigned int k = 0; k < opt->number_of_kernels; ++k) {
      vectorIterMult<<<grid, block>>>(d_vector, opt->number_of_threads,
                                      opt->inner_iterations);
    }
  }
#ifdef TIME_DETAILED
  time_elapsed[TOTAL_LAUNCH_COST] = getCpuSeconds() - time_start;
#endif
  cudaDeviceSynchronize();
#ifdef TIME_DETAILED
  time_elapsed[EXEC_TIME] = getCpuSeconds() - time_start;
#endif
  return 0;
}

static cudaGraph_t g_main_graph;
static cudaGraphNode_t *g_nodes;
static cudaGraphExec_t g_exec_work_graph;
static cudaStream_t g_stream_for_cuda_graph;

cudaError_t device_graph_setup(const struct options *opt, DataType **d_vector) {
#ifdef TIME_DETAILED
  double time_start = getCpuSeconds();
#endif
  cudaError_t device_error;
  device_error = cudaGraphCreate(&g_main_graph, 0);
  if (cudaSuccess != device_error) {
    printCudaError(device_error);
    return device_error;
  }

  dim3 block(TPB);
  dim3 grid((opt->number_of_threads + TPB - 1) / TPB);

  void *ka_kernel[] = {(void *)d_vector, (void *)&opt->number_of_threads,
                       (void *)&opt->inner_iterations};
  cudaKernelNodeParams np_kernel = {0};
  np_kernel.func = (void *)vectorIterMult;
  np_kernel.gridDim = grid;
  np_kernel.blockDim = block;
  np_kernel.kernelParams = ka_kernel;

  cudaGraphNode_t *last_node = NULL;
  unsigned int num_dependencies = 0;
  g_nodes = (cudaGraphNode_t *)malloc(opt->number_of_kernels *
                                      sizeof(cudaGraphNode_t));
  for (unsigned int i = 0; i < opt->number_of_kernels; ++i) {
    device_error = cudaGraphAddKernelNode(&g_nodes[i], g_main_graph, last_node,
                                          num_dependencies, &np_kernel);

    if (cudaSuccess != device_error) {
      printCudaError(device_error);
      return device_error;
    }
    last_node = &g_nodes[i];
    num_dependencies = 1;
  }

  device_error = cudaGraphInstantiateWithFlags(&g_exec_work_graph, g_main_graph,
                                      0);
  if (cudaSuccess != device_error) {
    printCudaError(device_error);
    return device_error;
  }

  device_error = cudaStreamCreateWithFlags(&g_stream_for_cuda_graph,
                                           cudaStreamNonBlocking);
  if (cudaSuccess != device_error) {
    printCudaError(device_error);
    return device_error;
  }
  device_error = cudaGraphUpload(g_exec_work_graph, g_stream_for_cuda_graph);
  if (cudaSuccess != device_error) {
    printCudaError(device_error);
    return device_error;
  }
#ifdef TIME_DETAILED
  time_elapsed[GRAPH_CREATION] = getCpuSeconds() - time_start;
#endif
  return cudaSuccess;
}

cudaError_t device_graph_run(const struct options *opt) {
  cudaError_t device_error;
#ifdef TIME_DETAILED
  double time_start;
  time_start = getCpuSeconds();
#endif
  for (unsigned int i = 0; i < opt->outer_iterations; ++i) {
    cudaGraphLaunch(g_exec_work_graph, g_stream_for_cuda_graph);
  }
#ifdef TIME_DETAILED
  time_elapsed[TOTAL_LAUNCH_COST] = getCpuSeconds() - time_start;
#endif
  cudaStreamSynchronize(g_stream_for_cuda_graph);
#ifdef TIME_DETAILED
  time_elapsed[EXEC_TIME] = getCpuSeconds() - time_start;
#endif
#ifdef MEM_CHECK
#if defined(__HIP)
  system("rocm-smi --showmeminfo vram");
#else
  system("nvidia-smi");
#endif
#endif
  return cudaSuccess;
}

void device_graph_teardown(void) {
  cudaStreamDestroy(g_stream_for_cuda_graph);
  cudaGraphExecDestroy(g_exec_work_graph);
  cudaGraphDestroy(g_main_graph);
  free(g_nodes);
}

