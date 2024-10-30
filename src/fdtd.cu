#include <stdbool.h>
#include <stdio.h>

#include "options.h"
#include "device.h"
#include "utils.h"
#include "host.h"
#include "em.h"

int run(const struct options *opt);

int main(int argc, char* argv[])
{
  struct options opts;
  int return_status;

  return_status = parse_arguments(&opts, argc, argv);
  if(0 != return_status){
    printf("Error parsing arguments.\n");
    return -1;
  }
  if (opts.print_options)
    print_options(&opts);

  if (opts.print_header){
    print_header();
    return 0;
  }

  return_status = run(&opts);

  return 0;
}

static double time_elapsed[TOTAL_NR_TIMES];

int run(const struct options *opt) {
  double time_start[2] = {0,0};
  reset_times(time_elapsed);
  time_start[0] = getCpuSeconds();

  /* Generate host data. */
  fields_t *h_fields;
  h_fields = host_setup(opt->Nx, opt->Ny, opt->Nz,
                        opt->Lx, opt->Ly, opt->Lz,
                        opt->seed, true);

  fields_t *d_fields;
  d_fields = device_setup(h_fields);

  if (opt->run_cpu) {
    printf("Running CPU version\n");
    cpu_kernel_run(h_fields, opt->timesteps, opt->sampling);
  }

  printf("Running %s version\n", opt->run_graph ? "Graph" : "Baseline");
  time_start[1] = getCpuSeconds();

  if (opt->run_graph) {
    device_graph_setup(opt, d_fields);
    device_graph_run(opt);
  } else {
    device_kernel_run(d_fields, opt->timesteps);
  }
  time_elapsed[CUDA_DIFF_TIME] = getCpuSeconds() - time_start[1];

  fields_t *r_fields;
  r_fields = host_setup(opt->Nx, opt->Ny, opt->Nz,
                        opt->Lx, opt->Ly, opt->Lz,
                        opt->seed, false);
  device_get_fields(r_fields,d_fields);

  if(opt->run_cpu){
    compaire_fields(h_fields,r_fields, 0.0000001);
  }

  /* teardown host resources. */
  host_teardown(h_fields);
  host_teardown(r_fields);

  /* teardown device resources. */
  device_teardown();
  if (opt->run_graph) {
    device_graph_teardown();
  }

  time_elapsed[TOTAL_TIME] = getCpuSeconds() - time_start[0];

  print_times(opt, time_elapsed, TOTAL_NR_TIMES);
  return 0;
}
