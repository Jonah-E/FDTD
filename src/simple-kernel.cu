#include "em.h"
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

  DataType *Hx = fields.Hx;
  DataType *Hy = fields.Hy;
  DataType *Hz = fields.Hz;

  DataType *Ex = fields.Ex;
  DataType *Ey = fields.Ey;
  DataType *Ez = fields.Ez;

  if (idx < Nx+1 && idy < Ny && idz < Nz) {
    Hx[idx + idy*(Nx+1) + idz*Ny*(Nx+1)] += (dt/MU0) * (
      (Ey[idx + idy*(Nx+1) + (idz+1)*Ny*(Nx+1)] -Ey[idx + idy*(Nx+1) + idz*Ny*(Nx+1)])/dz -
      (Ez[idx + (idy+1)*(Nx+1) + idz*(Ny+1)*(Nx+1)]-Ez[idx + idy*(Nx+1) + idz*(Ny+1)*(Nx+1)])/dy);
  }


  if (idx<Nx && idy < Ny+1 && idz < Nz){
    Hy[idx + idy*Nx + idz*(Ny+1)*Nx] += (dt/MU0) * (
      (Ez[(idx+1) + idy*(Nx+1) + idz*(Ny+1)*(Nx+1)]-Ez[idx + idy*(Nx+1) + idz*(Ny+1)*(Nx+1)])/dx -
      (Ex[idx + idy*Nx + (idz+1)*(Ny+1)*Nx]-Ex[idx + idy*Nx + idz*(Ny+1)*Nx])/dz);
  }

  if (idx<Nx && idy < Ny && idz < Nz+1){
    Hz[idx + idy*Nx + idz*Ny*Nx] += (dt/MU0) * (
      (Ex[idx + (idy+1)*Nx + idz*(Ny+1)*Nx]-Ex[idx + idy*Nx + idz*(Ny+1)*Nx])/dy -
      (Ey[(idx+1) + idy*(Nx+1) + idz*Ny*(Nx+1)]-Ey[idx + idy*(Nx+1) + idz*Ny*(Nx+1)])/dx);
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

  DataType *Hx = fields.Hx;
  DataType *Hy = fields.Hy;
  DataType *Hz = fields.Hz;

  DataType *Ex = fields.Ex;
  DataType *Ey = fields.Ey;
  DataType *Ez = fields.Ez;

  if (idx < Nx && idy < Ny && idz < Nz) {
    if (idy>0 && idz>0){
      Ex[idx + idy*Nx + idz*(Ny+1)*Nx] += (dt/EPS0) * (
          (Hz[idx+idy*Nx+idz*Ny*Nx]-Hz[idx+(idy-1)*Nx+idz*Ny*Nx])/dy -
          (Hy[idx+idy*Nx+idz*(Ny+1)*Nx]-Hy[idx+idy*Nx+(idz-1)*(Ny+1)*Nx])/dz);
    }
    if(idx>0 && idz>0){
      Ey[idx + idy*(Nx+1) + idz*Ny*(Nx+1)] += (dt/EPS0) * (
          (Hx[idx+idy*(Nx+1)+idz*Ny*(Nx+1)]-Hx[idx+idy*(Nx+1)+(idz-1)*Ny*(Nx+1)])/dz -
          (Hz[idx+idy*Nx+idz*Ny*Nx]-Hz[(idx-1)+idy*Nx+idz*Ny*Nx])/dx);
    }
    if(idx>0 && idy>0){
      Ez[idx + idy*(Nx+1) + idz*(Ny+1)*(Nx+1)] += (dt/EPS0) * (
          (Hy[idx+idy*Nx+idz*(Ny+1)*Nx]-Hy[(idx-1)+idy*Nx+idz*(Ny+1)*Nx])/dx -
          (Hx[idx+idy*(Nx+1)+idz*Ny*(Nx+1)]-Hx[idx+(idy-1)*(Nx+1)+idz*Ny*(Nx+1)])/dy);
    }
  }
}
