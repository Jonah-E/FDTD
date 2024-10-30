#include "options.h"
#include <argp.h>
#include <stdlib.h>

enum comp_opt {
  OPT_PRINT = 'p',
  OPT_PRINT_HEADER = 'h',
  OPT_GRAPH = 'g',
  OPT_CPU = 'c',
  OPT_TIMESTEP = 't',
  OPT_IT_BATCH = 'i',
  OPT_SAMPLING = 's',
  OPT_NX = 0x100,
  OPT_NY,
  OPT_NZ,
  OPT_LX,
  OPT_LY,
  OPT_LZ,
};

const char* argp_program_version = BUILD_VERSION;
static char doc[] = "";
static char args_doc[] = "";
static struct argp_option arguments[] = {
    {NULL,    OPT_PRINT, 0, 0, "Print options."},
    {NULL,    OPT_PRINT_HEADER, 0, 0, "Print csv header."},
    {"graph", OPT_GRAPH, 0, 0, "Execute graph version."},
    {"cpu",   OPT_CPU  , 0, 0, "Execute the reference cpu version (and compare.)"},
    {"timesteps", OPT_TIMESTEP, "timesteps", 0,
      "The number of timesteps, have to be a multiple of it_batch_size."},
    {"it_batch", OPT_IT_BATCH, "it_batch_size", 0, "Set the size of the iteration batch."},
    {"sampling", OPT_SAMPLING, "sampling", 0, "How often to sample the reference (cpu) version."},
    {"Nx", OPT_NX, "Nx", 0, "Number of cells in the x direction."},
    {"Ny", OPT_NY, "Ny", 0, "Number of cells in the y direction."},
    {"Nz", OPT_NZ, "Nz", 0, "Number of cells in the z direction."},
    {"Lx", OPT_LX, "Lx", 0, "The size of the cavity in the x direction."},
    {"Ly", OPT_LY, "Ly", 0, "The size of the cavity in the y direction."},
    {"Lz", OPT_LZ, "Lz", 0, "The size of the cavity in the z direction."},
    {0}};

static void reset_options(struct options* opt)
{
  opt->Nx = 25;
  opt->Ny = 20;
  opt->Nz = 15;
  opt->Lx = 0.05;
  opt->Ly = 0.04;
  opt->Lz = 0.03;
  opt->timesteps = 8000;
  opt->it_batch_size = 100;
  opt->sampling = 0;
  opt->run_graph = false;
  opt->run_cpu = false;
  opt->print_options = false;
  opt->print_header = false;
  opt->seed = 1;
}

static error_t parse_opt(int key, char* arg, struct argp_state* state)
{
  struct options* opt = state->input;

  switch (key) {
  case OPT_PRINT:
    opt->print_options = true;
    break;
  case OPT_PRINT_HEADER:
    opt->print_header = true;
    break;
  case OPT_GRAPH:
    opt->run_graph = true;
    break;
  case OPT_CPU:
    opt->run_cpu = true;
    break;
  case OPT_TIMESTEP:
    opt->timesteps = atoi(arg);
    break;
  case OPT_IT_BATCH:
    opt->it_batch_size = atoi(arg);
    break;
  case OPT_SAMPLING:
    opt->run_cpu = true;
    opt->sampling = atoi(arg);
    break;
  case OPT_NX:
      opt->Nx = atoi(arg);
    break;
  case OPT_NY:
      opt->Ny = atoi(arg);
    break;
  case OPT_NZ:
      opt->Nz = atoi(arg);
    break;
  case OPT_LX:
      opt->Lx = atof(arg);
    break;
  case OPT_LY:
      opt->Ly = atof(arg);
    break;
  case OPT_LZ:
      opt->Lz = atof(arg);
    break;
  case ARGP_KEY_ARG:
    break;
  default:
    return ARGP_ERR_UNKNOWN;
  }
  return 0;
}

static struct argp argp = {arguments, parse_opt, args_doc, doc, 0, 0, 0};

void print_options(const struct options* opt)
{
  printf("Options: \n"
         "\trun_graph = %s\n"
         "\trun_cpu = %s\n"
         "\timesteps = %d\n"
         "\tit_batch_size = %d\n"
         "\tNx = %d\n"
         "\tNy = %d\n"
         "\tNz = %d\n"
         "\tLx = %e\n"
         "\tLy = %e\n"
         "\tLz = %e\n",
          opt->run_graph ? "True" : "False",
          opt->run_cpu ? "True" : "False",
          opt->timesteps,
          opt->it_batch_size,
          opt->Nx,
          opt->Ny,
          opt->Nz,
          opt->Lx,
          opt->Ly,
          opt->Lz
         );
}

int parse_arguments(struct options* opt, int argc, char* argv[])
{
  error_t pars_error;
  reset_options(opt);
  pars_error = argp_parse(&argp, argc, argv, 0, 0, opt);

  if (pars_error)
    return -1;

  return 0;
}
