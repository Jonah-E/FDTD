#include "host.h"
#include "em.h"
#include "options.h"
#include "utils.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

static fields_t gh_fields;

fields_t *host_setup(int Nx, int Ny, int Nz,
                     DataType Lx, DataType Ly, DataType Lz, unsigned int seed) {
  gh_fields.Nx = Nx;
  gh_fields.Ny = Ny;
  gh_fields.Nz = Nz;
  gh_fields.Lx = Lx;
  gh_fields.Ly = Ly;
  gh_fields.Lz = Lz;
  gh_fields.dx = ((DataType) Lx)/Nx;
  gh_fields.dy = ((DataType) Ly)/Ny;
  gh_fields.dz = ((DataType) Lz)/Nz;

  gh_fields.dt = 1/(C0 * sqrt(
                      (1/(gh_fields.dx*gh_fields.dx)) +
                      (1/(gh_fields.dy*gh_fields.dy)) +
                      (1/(gh_fields.dz*gh_fields.dz))));

  gh_fields.Hx = (DataType *)calloc((Nx+1)*Ny*Nz, sizeof(DataType));
  gh_fields.Hy = (DataType *)calloc(Nx*(Ny+1)*Nz, sizeof(DataType));
  gh_fields.Hz = (DataType *)calloc(Nx*Ny*(Nz+1), sizeof(DataType));

  gh_fields.Ex = (DataType *)calloc(Nx*(Ny+1)*(Nz+1), sizeof(DataType));
  gh_fields.Ey = (DataType *)calloc((Nx+1)*Ny*(Nz+1), sizeof(DataType));
  gh_fields.Ez = (DataType *)calloc((Nx+1)*(Ny+1)*Nz, sizeof(DataType));

  /* Init E field randomly. */
  generateRandVector(gh_fields.Ex, Nx*(Ny+1)*(Nz+1), 0, 1.0, seed);
  generateRandVector(gh_fields.Ey, (Nx+1)*Ny*(Nz+1), 0, 1.0, seed);
  generateRandVector(gh_fields.Ez, (Nx+1)*(Ny+1)*Nz, 0, 1.0, seed);

  /*PEC boundry*/
  for(int x = 0; x < Nx+1; ++x){
    for(int y = 0; y < Ny+1; ++y){
      for(int z = 0; z < Nz+1; ++z){
        if (x < Nx && (y == 0 || z == 0 || y == Ny || z == Nz)){
          gh_fields.Ex[x + y*Nx  + z * Nx * Ny] = 0;
        }
        if ( y < Ny && (x == 0 || z == 0 || x == Nx || z == Nz)){
          gh_fields.Ey[x  + y*Nx + z * Nx * Ny] = 0;
        }
        if (z < Nz && ( x == 0 || y == 0 || x == Nx || y == Ny)){
          gh_fields.Ez[x  + y*Nx  + z * Nx * Ny] = 0;
        }
      }
    }
  }

  return &gh_fields;
}

void host_teardown(void) {
  free(gh_fields.Hx);
  free(gh_fields.Hy);
  free(gh_fields.Hz);

  free(gh_fields.Ex);
  free(gh_fields.Ey);
  free(gh_fields.Ez);
}

void update_fields(fields_t *h_fields) {

  /*Extract values for clearaty*/
  const int Nx = h_fields->Nx;
  const int Ny = h_fields->Ny;
  const int Nz = h_fields->Nz;

  const DataType dx = h_fields->dx;
  const DataType dy = h_fields->dy;
  const DataType dz = h_fields->dz;

  const DataType dt = h_fields->dt;

  DataType *Hx = h_fields->Hx;
  DataType *Hy = h_fields->Hy;
  DataType *Hz = h_fields->Hz;

  DataType *Ex = h_fields->Ex;
  DataType *Ey = h_fields->Ey;
  DataType *Ez = h_fields->Ez;

  /* Update the magnetic field*/
  for(int x=0; x < Nx+1; ++x){
    for(int y=0; y < Ny+1; ++y){
      for(int z=0; z < Nz+1; ++z){
        if (y < Ny && z < Nz) {
          Hx[x + y*Nx + z*Ny*Nx] += (dt/MU0) * (
            (Ey[x + y*Nx + (z+1)*Ny*Nx]-Ey[x + y*Nx + z*Ny*Nx])/dz -
            (Ez[x + (y+1)*Nx + z*Ny*Nx]-Ez[x + y*Nx + z*Ny*Nx])/dy);
        }

        if (x<Nx && z < Nz){
          Hy[x + y*Nx + z*Ny*Nx] += (dt/MU0) * (
            (Ez[(x+1) + y*Nx + z*Ny*Nx]-Ez[x + y*Nx + z*Ny*Nx])/dx -
            (Ex[x + y*Nx + (z+1)*Ny*Nx]-Ex[x + y*Nx + z*Ny*Nx])/dz);
        }

        if (x<Nx && y < Ny){
          Hz[x + y*Nx + z*Ny*Nx] += (dt/MU0) * (
            (Ex[x + (y+1)*Nx + z*Ny*Nx]-Ex[x + y*Nx + z*Ny*Nx])/dy -
            (Ey[(x+1) + y*Nx + z*Ny*Nx]-Ey[x + y*Nx + z*Ny*Nx])/dx);
        }
      }
    }
  }


  /* Update the electric field*/
  for(int x=0; x < Nx; ++x){
    for(int y=0; y < Ny; ++y){
      for(int z=0; z < Nz; ++z){
        if (y>0 && z>0){
          Ex[x + y*Nx + z*Ny*Nx] += (dt/EPS0) * (
              (Hz[x+y*Nx+z*Ny*Nx]-Hz[x+(y-1)*Nx+z*Ny*Nx])/dy -
              (Hy[x+y*Nx+z*Ny*Nx]-Hy[x+y*Nx+(z-1)*Ny*Nx])/dz);
        }
        if(x>0 && z>0){
          Ey[x + y*Nx + z*Ny*Nx] += (dt/EPS0) * (
              (Hx[x+y*Nx+z*Ny*Nx]-Hx[x+y*Nx+(z-1)*Ny*Nx])/dz -
              (Hz[x+y*Nx+z*Ny*Nx]-Hz[(x-1)+y*Nx+z*Ny*Nx])/dx);
        }
        if(x>0 && y>0){
          Ez[x + y*Nx + z*Ny*Nx] += (dt/EPS0) * (
              (Hy[x+y*Nx+z*Ny*Nx]-Hy[(x-1)+y*Nx+z*Ny*Nx])/dx -
              (Hx[x+y*Nx+z*Ny*Nx]-Hx[x+(y-1)*Nx+z*Ny*Nx])/dy);
        }
      }
    }
  }
}

int cpu_kernel_run(fields_t *h_fields, int time_steps, int sample_freq) {
  #define NR_SAMPLES 2
  int sample_loc[NR_SAMPLES][3] = {
    {h_fields->Nx/5, h_fields->Ny/5, h_fields->Nz/5},
    {h_fields->Nx/2, h_fields->Ny/2, h_fields->Nz/2}
  };

  if (sample_freq > 0){
    printf("tag,timestep");
    for (int i = 0; i < NR_SAMPLES; ++i){
      printf(",Ex_%d,Ey_%d,Ez_%d",i,i,i);
    }
    printf("\n");
  }
  for (int i = 0; i < time_steps; ++i){
    update_fields(h_fields);

    if (sample_freq > 0){
      if (i%sample_freq == 0){
        printf("tag,%d",i);
        for (int i = 0; i < NR_SAMPLES; ++i){
          int x,y,z;
          x = sample_loc[i][0];
          y = sample_loc[i][1];
          z = sample_loc[i][2];
          printf(",%f,%f,%f",
                  h_fields->Ex[x+y*h_fields->Nx+z*h_fields->Ny*h_fields->Nx],
                  h_fields->Ey[x+y*h_fields->Nx+z*h_fields->Ny*h_fields->Nx],
                  h_fields->Ez[x+y*h_fields->Nx+z*h_fields->Ny*h_fields->Nx]);

        }
        printf("\n");
      }
    }
  }

  return 0;
}
