#include "host.h"
#include "em.h"
#include "options.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

fields_t* host_setup(int Nx, int Ny, int Nz, DataType Lx, DataType Ly,
                     DataType Lz, unsigned int seed, bool generate)
{
  fields_t* gh_fields;

  gh_fields = (fields_t*) malloc(sizeof(fields_t));
  gh_fields->Nx = Nx;
  gh_fields->Ny = Ny;
  gh_fields->Nz = Nz;
  gh_fields->Lx = Lx;
  gh_fields->Ly = Ly;
  gh_fields->Lz = Lz;
  gh_fields->dx = ((DataType) Lx) / Nx;
  gh_fields->dy = ((DataType) Ly) / Ny;
  gh_fields->dz = ((DataType) Lz) / Nz;

  gh_fields->dt = 1 / (C0 * sqrt((1 / (gh_fields->dx * gh_fields->dx)) +
                                 (1 / (gh_fields->dy * gh_fields->dy)) +
                                 (1 / (gh_fields->dz * gh_fields->dz))));

  gh_fields->Hx = (DataType*) calloc((Nx + 1) * Ny * Nz, sizeof(DataType));
  gh_fields->Hy = (DataType*) calloc(Nx * (Ny + 1) * Nz, sizeof(DataType));
  gh_fields->Hz = (DataType*) calloc(Nx * Ny * (Nz + 1), sizeof(DataType));

  gh_fields->Ex =
      (DataType*) calloc(Nx * (Ny + 1) * (Nz + 1), sizeof(DataType));
  gh_fields->Ey =
      (DataType*) calloc((Nx + 1) * Ny * (Nz + 1), sizeof(DataType));
  gh_fields->Ez =
      (DataType*) calloc((Nx + 1) * (Ny + 1) * Nz, sizeof(DataType));

  if (generate) {

    /* Init E field randomly. */
    generateRandVector(gh_fields->Ex, Nx * (Ny + 1) * (Nz + 1), 0, 1.0, seed);
    generateRandVector(gh_fields->Ey, (Nx + 1) * Ny * (Nz + 1), 0, 1.0, seed);
    generateRandVector(gh_fields->Ez, (Nx + 1) * (Ny + 1) * Nz, 0, 1.0, seed);

    /*PEC boundry*/
    for (int x = 0; x < Nx + 1; ++x) {
      for (int y = 0; y < Ny + 1; ++y) {
        for (int z = 0; z < Nz + 1; ++z) {
          if (x < Nx && (y == 0 || z == 0 || y == Ny || z == Nz)) {
            gh_fields->Ex[x + y * Nx + z * Nx * (Ny + 1)] = 0;
          }
          if (y < Ny && (x == 0 || z == 0 || x == Nx || z == Nz)) {
            gh_fields->Ey[x + y * (Nx + 1) + z * (Nx + 1) * Ny] = 0;
          }
          if (z < Nz && (x == 0 || y == 0 || x == Nx || y == Ny)) {
            gh_fields->Ez[x + y * (Nx + 1) + z * (Nx + 1) * (Ny + 1)] = 0;
          }
        }
      }
    }
  }

  return gh_fields;
}

void host_teardown(fields_t* h_fields)
{
  free(h_fields->Hx);
  free(h_fields->Hy);
  free(h_fields->Hz);

  free(h_fields->Ex);
  free(h_fields->Ey);
  free(h_fields->Ez);

  free(h_fields);
}

static bool is_equal_f(DataType a, DataType b, DataType tolerance)
{
  return fabs(a - b) < tolerance;
}

bool compaire_fields(fields_t* field_1, fields_t* field_2, DataType tolerance)
{

  if (field_1->Nx != field_2->Nx) {
    printf("Field size Nx is not equal\n");
    return false;
  }
  if (field_1->Ny != field_2->Ny) {
    printf("Field size Ny is not equal\n");
    return false;
  }
  if (field_1->Nz != field_2->Nz) {
    printf("Field size Nz is not equal\n");
    return false;
  }
  if (!is_equal_f(field_1->Lx, field_2->Lx, tolerance)) {
    printf("Field size Lx is not equal\n");
    return false;
  }
  if (!is_equal_f(field_1->Ly, field_2->Ly, tolerance)) {
    printf("Field size Ly is not equal\n");
    return false;
  }
  if (!is_equal_f(field_1->Lz, field_2->Lz, tolerance)) {
    printf("Field size Lz is not equal\n");
    return false;
  }
  if (!is_equal_f(field_1->dx, field_2->dx, tolerance)) {
    printf("Field size dx is not equal\n");
    return false;
  }
  if (!is_equal_f(field_1->dy, field_2->dy, tolerance)) {
    printf("Field size dy is not equal\n");
    return false;
  }
  if (!is_equal_f(field_1->dz, field_2->dz, tolerance)) {
    printf("Field size dz is not equal\n");
    return false;
  }
  if (!is_equal_f(field_1->dt, field_2->dt, tolerance)) {
    printf("Field size dt is not equal\n");
    return false;
  }

  /*Extract values for clarity*/
  const int Nx = field_1->Nx;
  const int Ny = field_1->Ny;
  const int Nz = field_1->Nz;

  /* Compare h-fields*/
  for (int x = 0; x < Nx + 1; ++x) {
    for (int y = 0; y < Ny + 1; ++y) {
      for (int z = 0; z < Nz + 1; ++z) {
        if (y < Ny && z < Nz) {
          if (!is_equal_f(field_1->Hx[x + y * (Nx + 1) + z * Ny * (Nx + 1)],
                          field_2->Hx[x + y * (Nx + 1) + z * Ny * (Nx + 1)],
                          tolerance)) {
            printf("Field Hx is not equal at position (%d,%d,%d)\n", x, y, z);
            printf("\t %e != %e\n",
                   field_1->Hx[x + y * (Nx + 1) + z * Ny * (Nx + 1)],
                   field_2->Hx[x + y * (Nx + 1) + z * Ny * (Nx + 1)]);
            return false;
          }
        }

        if (x < Nx && z < Nz) {
          if (!is_equal_f(field_1->Hy[x + y * Nx + z * (Ny + 1) * Nx],
                          field_2->Hy[x + y * Nx + z * (Ny + 1) * Nx],
                          tolerance)) {
            printf("Field Hy is not equal at position (%d,%d,%d)\n", x, y, z);
            printf("\t %e != %e\n", field_1->Hy[x + y * Nx + z * (Ny + 1) * Nx],
                   field_2->Hy[x + y * Nx + z * (Ny + 1) * Nx]);
            return false;
          }
        }

        if (x < Nx && y < Ny) {
          if (!is_equal_f(field_1->Hz[x + y * Nx + z * Ny * Nx],
                          field_2->Hz[x + y * Nx + z * Ny * Nx], tolerance)) {
            printf("Field Hz is not equal at position (%d,%d,%d)\n", x, y, z);
            printf("\t %e != %e\n", field_1->Hz[x + y * Nx + z * Ny * Nx],
                   field_2->Hz[x + y * Nx + z * Ny * Nx]);
            return false;
          }
        }
      }
    }
  }

  /* Update the electric field*/
  for (int x = 0; x < Nx + 1; ++x) {
    for (int y = 0; y < Ny + 1; ++y) {
      for (int z = 0; z < Nz + 1; ++z) {
        if (x < Nx) {
          if (!is_equal_f(field_1->Ex[x + y * Nx + z * (Ny + 1) * Nx],
                          field_2->Ex[x + y * Nx + z * (Ny + 1) * Nx],
                          tolerance)) {
            printf("Field Ex is not equal at position (%d,%d,%d)\n", x, y, z);
            printf("\t %e != %e\n", field_1->Ex[x + y * Nx + z * (Ny + 1) * Nx],
                   field_2->Ex[x + y * Nx + z * (Ny + 1) * Nx]);
            return false;
          }
        }
        if (y < Ny) {
          if (!is_equal_f(field_1->Ey[x + y * (Nx + 1) + z * Ny * (Nx + 1)],
                          field_2->Ey[x + y * (Nx + 1) + z * Ny * (Nx + 1)],
                          tolerance)) {
            printf("Field Ey is not equal at position (%d,%d,%d)\n", x, y, z);
            printf("\t %e != %e\n",
                   field_1->Ey[x + y * (Nx + 1) + z * Ny * (Nx + 1)],
                   field_2->Ey[x + y * (Nx + 1) + z * Ny * (Nx + 1)]);
            return false;
          }
        }
        if (z < Nz) {
          if (!is_equal_f(
                  field_1->Ez[x + y * (Nx + 1) + z * (Ny + 1) * (Nx + 1)],
                  field_2->Ez[x + y * (Nx + 1) + z * (Ny + 1) * (Nx + 1)],
                  tolerance)) {
            printf("Field Ez is not equal at position (%d,%d,%d)\n", x, y, z);
            printf("\t %e != %e\n",
                   field_1->Ez[x + y * (Nx + 1) + z * (Ny + 1) * (Nx + 1)],
                   field_2->Ez[x + y * (Nx + 1) + z * (Ny + 1) * (Nx + 1)]);
            return false;
          }
        }
      }
    }
  }

  return true;
}

void update_fields(fields_t* h_fields)
{

  /*Extract values for clarity*/
  const int Nx = h_fields->Nx;
  const int Ny = h_fields->Ny;
  const int Nz = h_fields->Nz;

  const DataType dx = h_fields->dx;
  const DataType dy = h_fields->dy;
  const DataType dz = h_fields->dz;

  const DataType dt = h_fields->dt;

  DataType* Hx = h_fields->Hx;
  DataType* Hy = h_fields->Hy;
  DataType* Hz = h_fields->Hz;

  DataType* Ex = h_fields->Ex;
  DataType* Ey = h_fields->Ey;
  DataType* Ez = h_fields->Ez;

  /* Update the magnetic field*/
  for (int x = 0; x < Nx + 1; ++x) {
    for (int y = 0; y < Ny + 1; ++y) {
      for (int z = 0; z < Nz + 1; ++z) {
        if (y < Ny && z < Nz) {
          Hx[x + y * (Nx + 1) + z * Ny * (Nx + 1)] +=
              (dt / MU0) *
              ((Ey[x + y * (Nx + 1) + (z + 1) * Ny * (Nx + 1)] -
                Ey[x + y * (Nx + 1) + z * Ny * (Nx + 1)]) /
                   dz -
               (Ez[x + (y + 1) * (Nx + 1) + z * (Ny + 1) * (Nx + 1)] -
                Ez[x + y * (Nx + 1) + z * (Ny + 1) * (Nx + 1)]) /
                   dy);
        }

        if (x < Nx && z < Nz) {
          Hy[x + y * Nx + z * (Ny + 1) * Nx] +=
              (dt / MU0) *
              ((Ez[(x + 1) + y * (Nx + 1) + z * (Ny + 1) * (Nx + 1)] -
                Ez[x + y * (Nx + 1) + z * (Ny + 1) * (Nx + 1)]) /
                   dx -
               (Ex[x + y * Nx + (z + 1) * (Ny + 1) * Nx] -
                Ex[x + y * Nx + z * (Ny + 1) * Nx]) /
                   dz);
        }

        if (x < Nx && y < Ny) {
          Hz[x + y * Nx + z * Ny * Nx] +=
              (dt / MU0) * ((Ex[x + (y + 1) * Nx + z * (Ny + 1) * Nx] -
                             Ex[x + y * Nx + z * (Ny + 1) * Nx]) /
                                dy -
                            (Ey[(x + 1) + y * (Nx + 1) + z * Ny * (Nx + 1)] -
                             Ey[x + y * (Nx + 1) + z * Ny * (Nx + 1)]) /
                                dx);
        }
      }
    }
  }

  /* Update the electric field*/
  for (int x = 0; x < Nx; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {
        if (y > 0 && z > 0) {
          Ex[x + y * Nx + z * (Ny + 1) * Nx] +=
              (dt / EPS0) * ((Hz[x + y * Nx + z * Ny * Nx] -
                              Hz[x + (y - 1) * Nx + z * Ny * Nx]) /
                                 dy -
                             (Hy[x + y * Nx + z * (Ny + 1) * Nx] -
                              Hy[x + y * Nx + (z - 1) * (Ny + 1) * Nx]) /
                                 dz);
        }
        if (x > 0 && z > 0) {
          Ey[x + y * (Nx + 1) + z * Ny * (Nx + 1)] +=
              (dt / EPS0) * ((Hx[x + y * (Nx + 1) + z * Ny * (Nx + 1)] -
                              Hx[x + y * (Nx + 1) + (z - 1) * Ny * (Nx + 1)]) /
                                 dz -
                             (Hz[x + y * Nx + z * Ny * Nx] -
                              Hz[(x - 1) + y * Nx + z * Ny * Nx]) /
                                 dx);
        }
        if (x > 0 && y > 0) {
          Ez[x + y * (Nx + 1) + z * (Ny + 1) * (Nx + 1)] +=
              (dt / EPS0) * ((Hy[x + y * Nx + z * (Ny + 1) * Nx] -
                              Hy[(x - 1) + y * Nx + z * (Ny + 1) * Nx]) /
                                 dx -
                             (Hx[x + y * (Nx + 1) + z * Ny * (Nx + 1)] -
                              Hx[x + (y - 1) * (Nx + 1) + z * Ny * (Nx + 1)]) /
                                 dy);
        }
      }
    }
  }
}

int cpu_kernel_run(fields_t* h_fields, int time_steps, int sample_freq)
{
#define NR_SAMPLES 2
  int sample_loc[NR_SAMPLES][3] = {
      {h_fields->Nx / 5, h_fields->Ny / 5, h_fields->Nz / 5},
      {h_fields->Nx / 2, h_fields->Ny / 2, h_fields->Nz / 2}};

  if (sample_freq > 0) {
    printf("tag,timestep,dt");
    for (int i = 0; i < NR_SAMPLES; ++i) {
      printf(",Ex_%d,Ey_%d,Ez_%d", i, i, i);
    }
    printf("\n");
  }
  for (int i = 0; i < time_steps; ++i) {
    update_fields(h_fields);

    if (sample_freq > 0) {
      if (i % sample_freq == 0) {
        printf("tag,%d,%.10e", i, h_fields->dt);
        for (int i = 0; i < NR_SAMPLES; ++i) {
          int x, y, z;
          x = sample_loc[i][0];
          y = sample_loc[i][1];
          z = sample_loc[i][2];
          printf(",%e,%e,%e",
                 h_fields->Ex[x + y * h_fields->Nx +
                              z * (h_fields->Ny + 1) * h_fields->Nx],
                 h_fields->Ey[x + y * (h_fields->Nx + 1) +
                              z * h_fields->Ny * (h_fields->Nx + 1)],
                 h_fields->Ez[x + y * (h_fields->Nx + 1) +
                              z * (h_fields->Ny + 1) * (h_fields->Nx + 1)]);
        }
        printf("\n");
      }
    }
  }

  return 0;
}
