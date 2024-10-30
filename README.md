# FDTD

This is an 3D electromagnetic Finite-Difference Time-Domain (FDTD) simulation solving the 3D Maxwell equation for a vacuum cavity with a Perfect Electrical Conductor (PEC) boundary condition. The simulation is running on NVIDIA GPUs using CUDA with an option to either run a traditional implementation launching kernels in an loop, or an CUDA Graph version. The CUDA graph version bundles a number of kernel calls together in an Iteration Batch, which is then unrolled into a CUDA Graph. This enables the runtime to perform optimizations one the kernel scheduling and can speed-up the simulation.

# Build

To build the applicaiton you need to install the CUDA toolkit, [https://developer.nvidia.com/cuda-toolkit](https://developer.nvidia.com/cuda-toolkit). You also need gcc and make installed.

You can then build using the command:
````
make
````

# Running

The applicaiton can now be executed using

````
./fdtd
````

For option you can run

````
./fdtd --help
````


