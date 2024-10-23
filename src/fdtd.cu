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

  return_status = run(&opts);

  return 0;
}

static double time_elapsed[TOTAL_NR_TIMES];

int run(const struct options *opt) {
  cudaError_t device_error;
  double time_start[2] = {0,0};
  reset_times(time_elapsed);
  time_start[0] = getCpuSeconds();

  /* Generate host data. */
  fields_t *h_fields;
  h_fields = host_setup(opt->Nx, opt->Ny, opt->Nz,
                        opt->Lx, opt->Ly, opt->Lz,
                        opt->seed);
  if (NULL == h_fields) {
    return -1;
  }

  DataType *d_vector;

  if (opt->run_cpu) {
    cpu_kernel_run(h_fields, opt->timesteps, 1);
  } else {
    /* Setup device resources. */
    /*d_vector = device_setup(h_vector, opt->number_of_threads);
    if (NULL == d_vector) {
      device_teardown();
      host_teardown();
      return -1;
    }

    time_start[1] = getCpuSeconds();
    if (opt->run_graph) {
      device_error = device_graph_setup(opt, &d_vector);
      if (cudaSuccess == device_error) {
        device_graph_run(opt);
      }
    } else {
      device_kernel_run(opt, d_vector);
    }
    time_elapsed[CUDA_DIFF_TIME] = getCpuSeconds() - time_start[1];

    device_error = cudaMemcpy(h_vector, d_vector,
                              sizeof(DataType) * opt->number_of_threads,
                              cudaMemcpyDeviceToHost);
    if (device_error != cudaSuccess) {
      printCudaError(device_error);
    }*/

    /* teardown device resources. */
    /*
    device_teardown();
    if (opt->run_graph) {
      device_graph_teardown();
    }
    */
  }

  //DataType result = euclicianNormVector(h_vector, opt->number_of_threads);

  /* teardown host resources. */
  host_teardown();

  time_elapsed[TOTAL_TIME] = getCpuSeconds() - time_start[0];

  print_times(opt, time_elapsed, TOTAL_NR_TIMES, 0);
  return 0;
}
