# Toolchain file for Fugaku (A64FX, Fujitsu Clang via mpiFCCpx)
# Based on Makefile_target.inc Fugaku_CLANG section.

set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_SYSTEM_PROCESSOR aarch64)

set(CMAKE_CXX_COMPILER "mpiFCCpx" CACHE STRING "Fugaku MPI C++ compiler" FORCE)

# -Nclang: use Clang front-end of Fujitsu compiler
# -DFUJITSU_FX: macro expected by Bridge++ on Fugaku/A64FX
set(CMAKE_CXX_FLAGS_INIT "-Nclang -DFUJITSU_FX")
set(CMAKE_EXE_LINKER_FLAGS_INIT "-Nclang")

# High-optimisation flags (equivalent to Makefile OPT_HIGH)
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -mllvm -inline-threshold=1000"
    CACHE STRING "Release flags for Fugaku Clang" FORCE)
