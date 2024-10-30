#include "device.h"
#include "em.h"
#include <cassert>
#include <stdio.h>
#include <stdlib.h>

#define checkCudaErrors(cuda_returned_error_code)                              \
  {                                                                            \
    checkErrorCuda((cuda_returned_error_code), __FILE__, __LINE__);            \
  }

inline void checkErrorCuda(cudaError_t code, const char* file, int line)
{
  if (code) {
    fprintf(stderr, "CUDA Error: %s (%d) %s %d\n", cudaGetErrorString(code),
            code, file, line);
    exit(EXIT_FAILURE);
  }
}

__global__ void update_hfields(fields_t fields)
{
  const int idx = threadIdx.x + blockDim.x * blockIdx.x;
  const int idy = threadIdx.y + blockDim.y * blockIdx.y;
  const int idz = threadIdx.z + blockDim.z * blockIdx.z;

  /*Extract values for clarity*/
  const int Nx = fields.Nx;
  const int Ny = fields.Ny;
  const int Nz = fields.Nz;

  const DataType dx = fields.dx;
  const DataType dy = fields.dy;
  const DataType dz = fields.dz;

  const DataType dt = fields.dt;

  DataType* Hx = fields.Hx;
  DataType* Hy = fields.Hy;
  DataType* Hz = fields.Hz;

  DataType* Ex = fields.Ex;
  DataType* Ey = fields.Ey;
  DataType* Ez = fields.Ez;

  if (idx < Nx + 1 && idy < Ny && idz < Nz) {
    Hx[idx + idy * (Nx + 1) + idz * Ny * (Nx + 1)] +=
        (dt / MU0) *
        ((Ey[idx + idy * (Nx + 1) + (idz + 1) * Ny * (Nx + 1)] -
          Ey[idx + idy * (Nx + 1) + idz * Ny * (Nx + 1)]) /
             dz -
         (Ez[idx + (idy + 1) * (Nx + 1) + idz * (Ny + 1) * (Nx + 1)] -
          Ez[idx + idy * (Nx + 1) + idz * (Ny + 1) * (Nx + 1)]) /
             dy);
  }

  if (idx < Nx && idy < Ny + 1 && idz < Nz) {
    Hy[idx + idy * Nx + idz * (Ny + 1) * Nx] +=
        (dt / MU0) *
        ((Ez[(idx + 1) + idy * (Nx + 1) + idz * (Ny + 1) * (Nx + 1)] -
          Ez[idx + idy * (Nx + 1) + idz * (Ny + 1) * (Nx + 1)]) /
             dx -
         (Ex[idx + idy * Nx + (idz + 1) * (Ny + 1) * Nx] -
          Ex[idx + idy * Nx + idz * (Ny + 1) * Nx]) /
             dz);
  }

  if (idx < Nx && idy < Ny && idz < Nz + 1) {
    Hz[idx + idy * Nx + idz * Ny * Nx] +=
        (dt / MU0) * ((Ex[idx + (idy + 1) * Nx + idz * (Ny + 1) * Nx] -
                       Ex[idx + idy * Nx + idz * (Ny + 1) * Nx]) /
                          dy -
                      (Ey[(idx + 1) + idy * (Nx + 1) + idz * Ny * (Nx + 1)] -
                       Ey[idx + idy * (Nx + 1) + idz * Ny * (Nx + 1)]) /
                          dx);
  }
}

__global__ void update_efields(fields_t fields)
{
  const int idx = threadIdx.x + blockDim.x * blockIdx.x;
  const int idy = threadIdx.y + blockDim.y * blockIdx.y;
  const int idz = threadIdx.z + blockDim.z * blockIdx.z;

  /*Extract values for clearaty*/
  const int Nx = fields.Nx;
  const int Ny = fields.Ny;
  const int Nz = fields.Nz;

  const DataType dx = fields.dx;
  const DataType dy = fields.dy;
  const DataType dz = fields.dz;

  const DataType dt = fields.dt;

  DataType* Hx = fields.Hx;
  DataType* Hy = fields.Hy;
  DataType* Hz = fields.Hz;

  DataType* Ex = fields.Ex;
  DataType* Ey = fields.Ey;
  DataType* Ez = fields.Ez;

  if (idx < Nx && idy < Ny && idz < Nz) {
    if (idy > 0 && idz > 0) {
      Ex[idx + idy * Nx + idz * (Ny + 1) * Nx] +=
          (dt / EPS0) * ((Hz[idx + idy * Nx + idz * Ny * Nx] -
                          Hz[idx + (idy - 1) * Nx + idz * Ny * Nx]) /
                             dy -
                         (Hy[idx + idy * Nx + idz * (Ny + 1) * Nx] -
                          Hy[idx + idy * Nx + (idz - 1) * (Ny + 1) * Nx]) /
                             dz);
    }
    if (idx > 0 && idz > 0) {
      Ey[idx + idy * (Nx + 1) + idz * Ny * (Nx + 1)] +=
          (dt / EPS0) *
          ((Hx[idx + idy * (Nx + 1) + idz * Ny * (Nx + 1)] -
            Hx[idx + idy * (Nx + 1) + (idz - 1) * Ny * (Nx + 1)]) /
               dz -
           (Hz[idx + idy * Nx + idz * Ny * Nx] -
            Hz[(idx - 1) + idy * Nx + idz * Ny * Nx]) /
               dx);
    }
    if (idx > 0 && idy > 0) {
      Ez[idx + idy * (Nx + 1) + idz * (Ny + 1) * (Nx + 1)] +=
          (dt / EPS0) *
          ((Hy[idx + idy * Nx + idz * (Ny + 1) * Nx] -
            Hy[(idx - 1) + idy * Nx + idz * (Ny + 1) * Nx]) /
               dx -
           (Hx[idx + idy * (Nx + 1) + idz * Ny * (Nx + 1)] -
            Hx[idx + (idy - 1) * (Nx + 1) + idz * Ny * (Nx + 1)]) /
               dy);
    }
  }
}

static fields_t gd_fields;
fields_t* device_setup(fields_t* h_fields)
{
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

  checkCudaErrors(
      cudaMalloc(&gd_fields.Hx, sizeof(DataType) * (Nx + 1) * Ny * Nz));
  checkCudaErrors(
      cudaMalloc(&gd_fields.Hy, sizeof(DataType) * Nx * (Ny + 1) * Nz));
  checkCudaErrors(
      cudaMalloc(&gd_fields.Hz, sizeof(DataType) * Nx * Ny * (Nz + 1)));

  checkCudaErrors(
      cudaMalloc(&gd_fields.Ex, sizeof(DataType) * Nx * (Ny + 1) * (Nz + 1)));
  checkCudaErrors(
      cudaMalloc(&gd_fields.Ey, sizeof(DataType) * (Nx + 1) * Ny * (Nz + 1)));
  checkCudaErrors(
      cudaMalloc(&gd_fields.Ez, sizeof(DataType) * (Nx + 1) * (Ny + 1) * Nz));

  checkCudaErrors(cudaMemcpy(gd_fields.Hx, h_fields->Hx,
                             sizeof(DataType) * (Nx + 1) * Ny * Nz,
                             cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(gd_fields.Hy, h_fields->Hy,
                             sizeof(DataType) * Nx * (Ny + 1) * Nz,
                             cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(gd_fields.Hz, h_fields->Hz,
                             sizeof(DataType) * Nx * Ny * (Nz + 1),
                             cudaMemcpyHostToDevice));

  checkCudaErrors(cudaMemcpy(gd_fields.Ex, h_fields->Ex,
                             sizeof(DataType) * Nx * (Ny + 1) * (Nz + 1),
                             cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(gd_fields.Ey, h_fields->Ey,
                             sizeof(DataType) * (Nx + 1) * Ny * (Nz + 1),
                             cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(gd_fields.Ez, h_fields->Ez,
                             sizeof(DataType) * (Nx + 1) * (Ny + 1) * Nz,
                             cudaMemcpyHostToDevice));
  return &gd_fields;
}

void device_teardown(void)
{
  cudaFree(gd_fields.Hx);
  cudaFree(gd_fields.Hy);
  cudaFree(gd_fields.Hz);

  cudaFree(gd_fields.Ex);
  cudaFree(gd_fields.Ey);
  cudaFree(gd_fields.Ez);
}

#define TPB_X 8 // 16
#define TPB_Y 8 // 16
#define TPB_Z 8 // 16
int device_kernel_run(fields_t* d_fields, int timesteps)
{

  /*Total threads needs to be 1 more than Nx/Ny/Nz*/
  dim3 block(TPB_X, TPB_Y, TPB_Z);
  dim3 grid(((d_fields->Nx + 1 + TPB_X - 1) / TPB_X),
            ((d_fields->Ny + 1 + TPB_Y - 1) / TPB_Y),
            ((d_fields->Nz + 1 + TPB_Z - 1) / TPB_Z));

  for (unsigned int t = 0; t < timesteps; ++t) {
    update_hfields<<<grid, block>>>(*d_fields);
    checkCudaErrors(cudaGetLastError());
    update_efields<<<grid, block>>>(*d_fields);
    checkCudaErrors(cudaGetLastError());
  }

  cudaDeviceSynchronize();
  return 0;
}

void device_get_fields(fields_t* h_fields, fields_t* d_fields)
{
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

  checkCudaErrors(cudaMemcpy(h_fields->Hx, d_fields->Hx,
                             sizeof(DataType) * (Nx + 1) * Ny * Nz,
                             cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(h_fields->Hy, d_fields->Hy,
                             sizeof(DataType) * Nx * (Ny + 1) * Nz,
                             cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(h_fields->Hz, d_fields->Hz,
                             sizeof(DataType) * Nx * Ny * (Nz + 1),
                             cudaMemcpyDeviceToHost));

  checkCudaErrors(cudaMemcpy(h_fields->Ex, d_fields->Ex,
                             sizeof(DataType) * Nx * (Ny + 1) * (Nz + 1),
                             cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(h_fields->Ey, d_fields->Ey,
                             sizeof(DataType) * (Nx + 1) * Ny * (Nz + 1),
                             cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(h_fields->Ez, d_fields->Ez,
                             sizeof(DataType) * (Nx + 1) * (Ny + 1) * Nz,
                             cudaMemcpyDeviceToHost));
}

static cudaGraph_t g_main_graph;
static cudaGraphExec_t g_exec_work_graph;
static cudaStream_t g_stream_for_cuda_graph;

cudaError_t device_graph_setup(const struct options* opt, fields_t* d_fields)
{

  assert(((opt->timesteps % opt->it_batch_size) == 0) &&
         "Timestep have to be a mulitple of iteration batch size");

  checkCudaErrors(cudaGraphCreate(&g_main_graph, 0));

  /*Total threads needs to be 1 more than Nx/Ny/Nz*/
  dim3 block(TPB_X, TPB_Y, TPB_Z);
  dim3 grid(((d_fields->Nx + 1 + TPB_X - 1) / TPB_X),
            ((d_fields->Ny + 1 + TPB_Y - 1) / TPB_Y),
            ((d_fields->Nz + 1 + TPB_Z - 1) / TPB_Z));

  /* H-fields kernel node*/
  void* ka_kernel_h[] = {(void*) d_fields};
  cudaKernelNodeParams np_kernel_h = {0};
  np_kernel_h.func = (void*) update_hfields;
  np_kernel_h.gridDim = grid;
  np_kernel_h.blockDim = block;
  np_kernel_h.kernelParams = ka_kernel_h;

  /* E-fields kernel node*/
  void* ka_kernel_e[] = {(void*) d_fields};
  cudaKernelNodeParams np_kernel_e = {0};
  np_kernel_e.func = (void*) update_efields;
  np_kernel_e.gridDim = grid;
  np_kernel_e.blockDim = block;
  np_kernel_e.kernelParams = ka_kernel_e;

  cudaGraphNode_t* last_node_p = NULL;
  cudaGraphNode_t last_node, current_node;
  unsigned int num_dependencies = 0;

  for (unsigned int i = 0; i < opt->it_batch_size; ++i) {
    checkCudaErrors(cudaGraphAddKernelNode(&current_node, g_main_graph,
                                           last_node_p, num_dependencies,
                                           &np_kernel_h));

    last_node = current_node;
    last_node_p = &last_node;
    num_dependencies = 1;

    checkCudaErrors(cudaGraphAddKernelNode(&current_node, g_main_graph,
                                           last_node_p, num_dependencies,
                                           &np_kernel_e));

    last_node = current_node;
    last_node_p = &last_node;
    num_dependencies = 1;
  }

  checkCudaErrors(
      cudaGraphInstantiateWithFlags(&g_exec_work_graph, g_main_graph, 0));

  checkCudaErrors(cudaStreamCreateWithFlags(&g_stream_for_cuda_graph,
                                            cudaStreamNonBlocking));
  checkCudaErrors(cudaGraphUpload(g_exec_work_graph, g_stream_for_cuda_graph));

  return cudaSuccess;
}

cudaError_t device_graph_run(const struct options* opt)
{
  // Timestep is checked in device_graph_setup(..)
  int nr_graph_launches = opt->timesteps / opt->it_batch_size;

  for (unsigned int i = 0; i < nr_graph_launches; ++i) {
    cudaGraphLaunch(g_exec_work_graph, g_stream_for_cuda_graph);
    checkCudaErrors(cudaGetLastError());
  }
  cudaStreamSynchronize(g_stream_for_cuda_graph);
  return cudaSuccess;
}

void device_graph_teardown(void)
{
  cudaStreamDestroy(g_stream_for_cuda_graph);
  cudaGraphExecDestroy(g_exec_work_graph);
  cudaGraphDestroy(g_main_graph);
}
