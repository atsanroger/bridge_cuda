cmake_minimum_required(VERSION 3.21)
project(Bridge VERSION 2.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ----------- Options (mirrors Makefile variables) -----------
set(USE_MPI "ON" CACHE STRING "MPI mode: ON | OFF | fjmpi")
set_property(CACHE USE_MPI PROPERTY STRINGS ON OFF fjmpi)
option(USE_OPENMP      "Enable OpenMP"                 ON)
option(USE_OPT_CODE    "Use optimised (Imp) variants"  ON)
option(USE_ALTERNATIVE "Enable alt-code (lib_alt)"     ON)
option(USE_ALT_ACCEL   "Enable accelerator code"       OFF)

set(USE_ACCEL_IMPL "cuda" CACHE STRING "Accelerator backend: cuda | openacc | nvshmem")
set_property(CACHE USE_ACCEL_IMPL PROPERTY STRINGS cuda openacc nvshmem)

option(USE_DIRECT_COMM "Enable direct GPU-to-GPU communication (CUDA only)" OFF)

option(USE_ALT_QXS  "Enable QXS alt-code (lib_alt_QXS, Fujitsu A64FX)" OFF)
set(USE_QXS_ARCH "general" CACHE STRING "QXS SIMD arch: acle | general")
set_property(CACHE USE_QXS_ARCH PROPERTY STRINGS acle general)

option(USE_LIME     "Enable LIME I/O library"              OFF)
option(USE_SFMT     "Enable SFMT random number library"    OFF)
option(USE_NTL      "Enable NTL number theory library"     OFF)
option(USE_FFTW     "Enable FFTW library"                  OFF)
option(USE_YAMLCPP  "Enable yaml-cpp parameter library"    OFF)
option(USE_TINYXML2 "Enable tinyxml2 XML library"          OFF)
option(USE_QWS      "Enable QWS library"                   OFF)

set(BRIDGE_GAUGE_GROUP "SU3" CACHE STRING "Gauge group: SU3 | SU2 | SU_N")
set_property(CACHE BRIDGE_GAUGE_GROUP PROPERTY STRINGS SU3 SU2 SU_N)

# ----------- CMake helper modules (cmake/) -----------
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Source-tree symlinks (BridgeACC backend + QXS inline arch), mirroring the
# Makefile build's link-inline target.  See cmake/BridgeSymlinks.cmake.
include(BridgeSymlinks)

# ----------- CUDA/NVSHMEM sub-library (must precede src/ so bridgecuda is defined first) -----------
if(USE_ALT_ACCEL AND (USE_ACCEL_IMPL STREQUAL "cuda" OR USE_ACCEL_IMPL STREQUAL "nvshmem"))
    enable_language(CUDA)
    if(USE_ACCEL_IMPL STREQUAL "cuda")
        add_subdirectory(src/lib_alt_Accel/BridgeCUDA)
    else()
        add_subdirectory(src/lib_alt_Accel/BridgeNVSHMEM)
    endif()
endif()

# ----------- Main library -----------
add_subdirectory(src)

add_subdirectory(tests)

# ----------- Standalone test_alt_* programs -----------
# (the foreach + the batched-AMG bridgecuda_mrhs link).  See cmake/BridgeAltTests.cmake.
include(BridgeAltTests)
