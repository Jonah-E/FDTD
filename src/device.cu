#include "device.h"
#include "em.h"
#include "host.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>

#define checkCudaErrors(cuda_returned_error_code)                               \
  { checkErrorCuda((cuda_returned_error_code), __FILE__, __LINE__); }

inline void checkErrorCuda(cudaError_t code, const char *file, int line) {
  if (code) {
    fprintf(stderr, "CUDA Error: %s (%d) %s %d\n", cudaGetErrorString(code), code,
            file, line);
    exit(EXIT_FAILURE);
  }
}

#include "simple-kernel.cu"

static fields_t gd_fields;
fields_t *device_setup(fields_t *h_fields) {
  gd_fields.Nx = h_fields->Nx;
  gd_fields.Ny = h_fields->Ny;
  gd_fields.Nz = h_fields->Nz;
  gd_fields.Lx = h_fields->Lx;
  gd_fields.Ly = h_fields->Ly;
  gd_fields.Lz = h_fields->Lz;
  gd_fields.dx = h_fields->dx;
  gd_fields.dy = h_fields->dy;
  gd_fields.dz = h_fields->dz;

  gd_fields.dt = h_fields->dt;

  /*Extract sizes for clarity.*/
  int Nx = gd_fields.Nx;
  int Ny = gd_fields.Ny;
  int Nz = gd_fields.Nz;

  checkCudaErrors(cudaMalloc(&gd_fields.Hx, sizeof(DataType) * (Nx+1)*Ny*Nz));
  checkCudaErrors(cudaMalloc(&gd_fields.Hy, sizeof(DataType) * Nx*(Ny+1)*Nz));
  checkCudaErrors(cudaMalloc(&gd_fields.Hz, sizeof(DataType) * Nx*Ny*(Nz+1)));

  checkCudaErrors(cudaMalloc(&gd_fields.Ex, sizeof(DataType) * Nx*(Ny+1)*(Nz+1)));
  checkCudaErrors(cudaMalloc(&gd_fields.Ey, sizeof(DataType) * (Nx+1)*Ny*(Nz+1)));
  checkCudaErrors(cudaMalloc(&gd_fields.Ez, sizeof(DataType) * (Nx+1)*(Ny+1)*Nz));

  checkCudaErrors(cudaMemcpy(gd_fields.Hx, h_fields->Hx, sizeof(DataType) * (Nx+1)*Ny*Nz,
                             cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(gd_fields.Hy, h_fields->Hy, sizeof(DataType) * Nx*(Ny+1)*Nz,
                             cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(gd_fields.Hz, h_fields->Hz, sizeof(DataType) * Nx*Ny*(Nz+1),
                             cudaMemcpyHostToDevice));

  checkCudaErrors(cudaMemcpy(gd_fields.Ex, h_fields->Ex, sizeof(DataType) * Nx*(Ny+1)*(Nz+1),
                             cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(gd_fields.Ey, h_fields->Ey, sizeof(DataType) * (Nx+1)*Ny*(Nz+1),
                             cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(gd_fields.Ez, h_fields->Ez, sizeof(DataType) * (Nx+1)*(Ny+1)*Nz,
                             cudaMemcpyHostToDevice));
  return &gd_fields;
}

void device_teardown(void) {
  cudaFree(gd_fields.Hx);
  cudaFree(gd_fields.Hy);
  cudaFree(gd_fields.Hz);

  cudaFree(gd_fields.Ex);
  cudaFree(gd_fields.Ey);
  cudaFree(gd_fields.Ez);
}


#define TPB_X 8//16
#define TPB_Y 8//16
#define TPB_Z 8//16
int device_kernel_run(fields_t *d_fields, int timesteps) {

  /*Total threads needs to be 1 more than Nx/Ny/Nz*/
  dim3 block(TPB_X,TPB_Y,TPB_Z);
  dim3 grid(
    ((d_fields->Nx + 1 + TPB_X - 1) / TPB_X),
    ((d_fields->Ny + 1 + TPB_Y - 1) / TPB_Y),
    ((d_fields->Nz + 1 + TPB_Z - 1) / TPB_Z));


#ifdef TIME_DETAILED
  double time_start;
  time_start = getCpuSeconds();
#endif
  for (unsigned int t = 0; t < timesteps; ++t) {
    update_hfields<<<grid,block>>>(*d_fields);
    checkCudaErrors(cudaGetLastError());
    update_efields<<<grid,block>>>(*d_fields);
    checkCudaErrors(cudaGetLastError());
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

void device_get_fields(fields_t *h_fields, fields_t *d_fields){
  /*Extract sizes for clarity.*/
  int Nx = d_fields->Nx;
  int Ny = d_fields->Ny;
  int Nz = d_fields->Nz;

  h_fields->Nx = d_fields->Nx;
  h_fields->Ny = d_fields->Ny;
  h_fields->Nz = d_fields->Nz;
  h_fields->Lx = d_fields->Lx;
  h_fields->Ly = d_fields->Ly;
  h_fields->Lz = d_fields->Lz;
  h_fields->dx = d_fields->dx;
  h_fields->dy = d_fields->dy;
  h_fields->dz = d_fields->dz;

  h_fields->dt = d_fields->dt;

  checkCudaErrors(cudaMemcpy(h_fields->Hx, d_fields->Hx, sizeof(DataType) * (Nx+1)*Ny*Nz,
                             cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(h_fields->Hy, d_fields->Hy, sizeof(DataType) * Nx*(Ny+1)*Nz,
                             cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(h_fields->Hz, d_fields->Hz, sizeof(DataType) * Nx*Ny*(Nz+1),
                             cudaMemcpyDeviceToHost));

  checkCudaErrors(cudaMemcpy(h_fields->Ex, d_fields->Ex, sizeof(DataType) * Nx*(Ny+1)*(Nz+1),
                             cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(h_fields->Ey, d_fields->Ey, sizeof(DataType) * (Nx+1)*Ny*(Nz+1),
                             cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(h_fields->Ez, d_fields->Ez, sizeof(DataType) * (Nx+1)*(Ny+1)*Nz,
                             cudaMemcpyDeviceToHost));
}

static cudaGraph_t g_main_graph;
static cudaGraphNode_t *g_nodes;
static cudaGraphExec_t g_exec_work_graph;
static cudaStream_t g_stream_for_cuda_graph;

cudaError_t device_graph_setup(const struct options *opt, DataType **d_vector) {
#ifdef TIME_DETAILED
  double time_start = getCpuSeconds();
#endif
  /*
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
*/
  return cudaSuccess;
}

cudaError_t device_graph_run(const struct options *opt) {
#ifdef TIME_DETAILED
  double time_start;
  time_start = getCpuSeconds();
#endif
  for (unsigned int i = 0; i < 10; ++i) {
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

