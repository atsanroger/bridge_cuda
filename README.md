# Bridge C++ Codebase CMake & AMG Usage Guide

## 1. How to Use CMake

In this project, using CMake to configure and build is straightforward. We rely on CMake Presets to simplify the build commands for different environments and hardware architectures.

### Configure
Use the following command to generate the build files (Makefile or Ninja) using a preset:
```bash
cmake --preset <preset_name>
```

### Build
After configuring, use the following command to compile the project:
```bash
cmake --build --preset <preset_name>
```
This automatically uses the parallel build jobs (e.g., `jobs: 8`) defined in the preset.

## 2. CMake Presets Options

The `CMakePresets.json` in the project defines build configurations for various environments. They are categorized into the following series:

### CUDA Series (For NVIDIA GPUs)
- `cuda-ampere-mpi`: For A100 / A6000 (sm_80/86), uses MPI + OpenMP (Release build, for production)
- `cuda-ampere-debug`: For A100 / A6000 (sm_80/86), uses MPI (Debug build)
- `cuda-hopper-mpi`: For H100 (sm_90), uses MPI + OpenMP (Release build, for production)

### NVSHMEM Series (Multi-GPU Communication Optimization)
- `nvshmem-ampere-mpi`: For A100 (sm_80), combining CUDA with NVSHMEM, uses MPI + OpenMP
- `nvshmem-hopper-mpi`: For H100 (sm_90), combining CUDA with NVSHMEM, uses MPI + OpenMP

### CPU (PC) Series
- `pc-gnu-mpi`: CPU-only computation using GCC, Release build (MPI + OpenMP)
- `pc-gnu-debug`: CPU-only computation using GCC, Debug build

### OpenACC Series (OpenACC GPU Implementation)
- `openacc-ampere-mpi`: For A100 (cc80), compiled with NVHPC `nvc++`, uses MPI + OpenMP
- `openacc-hopper-mpi`: For H100 (cc90), compiled with NVHPC `nvc++`, uses MPI + OpenMP

### Other Platforms
- `fugaku`: Configuration for the Fugaku supercomputer, uses Clang + FJMPI

## 3. How to Start AMG

> [!WARNING]
> **Important Note**: Currently, AMG (Algebraic Multigrid) is only supported on **CUDA**. Please ensure you are building with a CUDA preset (e.g., `cuda-ampere-mpi`) when running AMG tests.

The tests and executables for AMG are located in test directories such as `test_alt_amg_dw_GEMM` and `test_alt_dw_multigrid`.

The basic steps to start and run AMG are as follows:

1. **Build the Test**:
   Make sure you have already built the base library via CMake. Then go to an AMG test directory (e.g., `test_alt_amg_dw_GEMM`). If you need to compile the test separately, run:
   ```bash
   cd test_alt_amg_dw_GEMM
   make
   ```
   *(Note: This typically generates an executable named `bridge_test_alt` inside the `build/` directory.)*

2. **Enter the Run Directory**:
   Tests are usually executed inside the `run` folder, which contains the necessary configuration files:
   ```bash
   cd test_alt_amg_dw_GEMM/run
   ```

3. **Run the Executable**:
   Inside the `run` directory, there is usually an executable named `exe` or a symlink to the built `bridge_test_alt`. 
   Run it alongside a YAML configuration file via MPI (or directly, depending on your environment):
   ```bash
   # For single GPU / direct execution
   ./exe test_alt_Multigrid.yaml
   
   # For multi-GPU execution using MPI
   mpirun -np <number_of_processes> ./exe test_alt_Multigrid.yaml
   ```
   *Note: You can adjust the AMG grid levels, solver parameters, and other settings by editing the YAML file. The `run` directory also contains historical log files (e.g., `run_amg.log`) that you can reference.*
