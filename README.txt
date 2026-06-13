Release Note for Bridge++ ver.2.0.0                 01 Mar 2023
====================================

   *************************************************************
     This package is free software available under
     "the GNU General Public License" published by the FSF.
   *************************************************************

Overview:
---------

Bridge++ is a code set for lattice gauge theory simulations
on Linux/Unix workstations and supercomputers (Fugaku,
FX1000:Fujitsu, SX-Aurora TSUBASA:NEC, etc) using "C/C++"
standard language with MPI for communication and OpenMP for
multi-threading.


Software Compatibility:
-----------------------

This code set follows the standard specification of C++.
The code has been tested with GNU C/C++ 9.4, Clang C++ 14.0,
Intel C++ 2020.4 and 2023.0, NVIDIA HPC SDK C++ 22.7,
Fujitsu C++ (Clang mode) 1.2.36, NEC SX C++ 3.5.1.


Installation:
-------------

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


Getting Help:
-------------

Visit our website, http://bridge.kek.jp/Lattice-code/ for a
contact information. See also our first step guide and the
doxygen document.

Bug reports are welcome. Please send them to the contact
address on the above website.


Acknowledging Bridge++:
-----------------------

If you use this software for your research, please cite
J.Phys.Conf.Ser. 523 (2014) 012046
J.Phys.Conf.Ser. 2207 (2022) 012053
as well as our website.

For example, "This work is in part based on Bridge++ code
(http://bridge.kek.jp/Lattice-code/)~\cite{Ueda:2014rya,Akahoshi:2021gvk}".


Project Member:
---------------
Tatsumi   Aoyama    (Univ. Tokyo)
Issaku    Kanamori  (RIKEN R-CCS)
Kazuyuki  Kanaya    (Univ. of Tsukuba)
Hideo     Matsufuru (KEK)
Yusuke    Namekawa  (Hiroshima Univ.)
Hidekatsu Nemura    (Kyoto Univ.)
Yusuke    Taniguchi (Univ. of Tsukuba)

Contributed by:
---------------
Sinya     Aoki      (Kyoto Univ.)
Takumi    Doi       (RIKEN)
Shoji     Hashimoto (KEK)
Noriyoshi Ishii     (Osaka Univ.)
Ken-ichi  Ishikawa  (Hiroshima Univ.)
Takashi   Kaneko    (KEK)
Yoshinobu Kuramashi (Univ. of Tsukuba)
Keigo     Nitadori  (RIKEN R-CCS)
Kenji     Sasaki    (Osaka Univ.)
Naoya     Ukita     (Univ. of Tsukuba)
Tomoteru  Yoshie    (Univ. of Tsukuba)

Former Member:
---------------
Yutaro    Akahoshi
Guido     Cossu
Takaya    Miyamoto
Shinji    Motoki
Jun-Ichi  Noaki
Kenji     Ogawa
Hana      Saito
Satoru    Ueda


This work is supported by the following grants, research programs,
and organizations:
- Grant-in-Aid for Scientific Research on Innovative Areas (of
  the Japanese Ministry of Education, Culture, Sports, Science
  and Technology) (Nos.20105001, 20105005)
- JSPS KAKENHI Grant Numbers 25400284, 15K05068, 16K05340,
  16H03988, 20K03961.
- HPCI Strategic Program Field 5 ("The origin of matter and
  the universe"),
- Priority Issue on Post-K computer (Elucidation of the
  Fundamental Laws and Evolution of the Universe)
- Joint Institute for Computational Fundamental Science (JICFuS)
- Program for Promoting Researches on the Supercomputer Fugaku
  "Simulation for basic science: from fundamental laws of
   particles to creation of nuclei"
- Multidisciplinary Cooperative Research Program, Center for
  Computational Science, University of Tsukuba (15a33, 16a42,
  17a41, 18-14, 19-43, 20-79, 21-42, 22-60)
- KEK Large-scale Simulation Program (Ohgata 09-22, 09/10-23,
  10-18, (T)11-09, 12-18, 12/13-15)
- Particle, Nuclear, and Astro Physics Simulation Program,
  Institute of Particle and Nuclear Studies, KEK (2019-T003,
  2019-004, 2020-002, 2021-006, 2022-002)
- KEK Computing Research Center
- Research Center for Nuclear Physics, Osaka University
- Yukawa Institute for Theoretical Physics, Kyoto University
