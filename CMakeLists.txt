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

# ----------- BridgeACC symlink (mirrors Makefile link-inline target) -----------
# Source files use #include "lib_alt_Accel/BridgeACC/..." which requires this symlink.
if(USE_ALT_ACCEL)
    if(USE_ACCEL_IMPL STREQUAL "cuda")
        set(_BRIDGE_ACC_TARGET "BridgeCUDA")
    elseif(USE_ACCEL_IMPL STREQUAL "nvshmem")
        set(_BRIDGE_ACC_TARGET "BridgeNVSHMEM")
    elseif(USE_ACCEL_IMPL STREQUAL "openacc")
        set(_BRIDGE_ACC_TARGET "BridgeOpenACC")
    endif()
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E create_symlink "${_BRIDGE_ACC_TARGET}" "BridgeACC"
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/src/lib_alt_Accel"
    )
endif()

# ----------- QXS inline symlink (mirrors Makefile link-inline for lib_alt_QXS) -----------
# lib_alt_QXS sources use #include "lib_alt_QXS/inline/..." which requires this symlink.
if(USE_ALT_QXS)
    if(USE_QXS_ARCH STREQUAL "acle")
        set(_QXS_INLINE "inline_ACLE")
    else()
        set(_QXS_INLINE "inline_General")
    endif()
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E create_symlink "${_QXS_INLINE}" "inline"
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/src/lib_alt_QXS"
    )
endif()

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

# ----------- test_alt_* programs -----------
file(READ "${CMAKE_SOURCE_DIR}/VERSION" _BRIDGE_VERSION)
string(STRIP "${_BRIDGE_VERSION}" _BRIDGE_VERSION)

foreach(_dir IN ITEMS
        test_alt_multigrid
        test_alt_dw_multigrid
        test_alt_dw_multigrid_dev2
        test_alt_dw_multigrid_dev4
        test_alt_dw_multigrid_dev5
        test_alt_dw_multigrid_dev6
        test_alt_spectrum
        test_alt_smear
        test_alt_hmc)
    set(_root_dir "${CMAKE_SOURCE_DIR}/${_dir}")
    if(NOT IS_DIRECTORY "${_root_dir}")
        continue()
    endif()

    # Flat layout: *.cpp directly in root (test_alt_multigrid, test_alt_dw_multigrid*)
    file(GLOB _srcs CONFIGURE_DEPENDS "${_root_dir}/*.cpp")
    if(_srcs)
        set(_inc_dir "${_root_dir}")
    else()
        # src/ layout: test_alt_spectrum, test_alt_smear, test_alt_hmc
        file(GLOB _srcs CONFIGURE_DEPENDS
            "${_root_dir}/src/*.cpp"
            "${_root_dir}/src/Tests/*.cpp")
        set(_inc_dir "${_root_dir}/src")
    endif()

    if(NOT _srcs)
        continue()
    endif()

    string(REPLACE "test_alt_" "" _short "${_dir}")
    add_executable(bridge_test_alt_${_short} ${_srcs})
    target_include_directories(bridge_test_alt_${_short} PRIVATE "${_inc_dir}")
    target_compile_definitions(bridge_test_alt_${_short} PRIVATE
        "_BRIDGE_VERSION=\"${_BRIDGE_VERSION}\"")
    target_link_libraries(bridge_test_alt_${_short} PRIVATE bridge)
endforeach()
