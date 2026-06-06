# ============================================================
#  Bridge main library  →  libbridge.a
# ============================================================

# ----------- variant selection -----------
if(USE_OPT_CODE)
    set(_OPT "Imp")
else()
    set(_OPT "Org")
endif()

if(USE_MPI STREQUAL "fjmpi")
    set(_COMMUN "FJMPI")
elseif(USE_MPI)
    set(_COMMUN "MPI")
else()
    set(_COMMUN "Single")
endif()

if(USE_OPENMP)
    set(_THREAD "OpenMP")
else()
    set(_THREAD "OpenMP_stub")
endif()

# ----------- source collection (CONFIGURE_DEPENDS: re-glob on new files) -----------

# Core lib: both Fopr/Org and Fopr/Imp are always included (factory pattern uses both)
file(GLOB BRIDGE_CORE CONFIGURE_DEPENDS
    "lib/*.cpp"
    "lib/Action/*.cpp"
    "lib/Action/Fermion/*.cpp"
    "lib/Action/Gauge/*.cpp"
    "lib/Eigen/*.cpp"
    "lib/Field/*.cpp"
    "lib/Fopr/*.cpp"
    "lib/Fopr/Org/*.cpp"
    "lib/Fopr/Imp/*.cpp"
    "lib/Force/Fermion/*.cpp"
    "lib/Force/Gauge/*.cpp"
    "lib/HMC/*.cpp"
    "lib/IO/*.cpp"
    "lib/Measurements/Fermion/*.cpp"
    "lib/Measurements/Gauge/*.cpp"
    "lib/Parameters/*.cpp"
    "lib/Smear/*.cpp"
    "lib/Solver/*.cpp"
    "lib/Tools/*.cpp"
)

# Variant dirs (only one of each pair compiled)
file(GLOB BRIDGE_COMMUN CONFIGURE_DEPENDS "lib/Communicator/${_COMMUN}/*.cpp")
file(GLOB BRIDGE_OPT    CONFIGURE_DEPENDS
    "lib/Field/${_OPT}/*.cpp"
    "lib/Force/Fermion/${_OPT}/*.cpp"
    "lib/Measurements/Fermion/${_OPT}/*.cpp"
)
file(GLOB BRIDGE_THREAD CONFIGURE_DEPENDS "lib/ResourceManager/${_THREAD}/*.cpp")

# lib_alt: entirely header-only templates except bridge_alt_init.cpp
set(BRIDGE_ALT_SOURCES)
if(USE_ALTERNATIVE)
    file(GLOB BRIDGE_ALT_SOURCES CONFIGURE_DEPENDS "lib_alt/*.cpp")
endif()

# lib_alt_Accel: CPU-side C++ wrappers (not the CUDA kernels, those are in libbridgecuda.a)
set(BRIDGE_ACCEL_SOURCES)
if(USE_ALT_ACCEL)
    file(GLOB BRIDGE_ACCEL_SOURCES CONFIGURE_DEPENDS
        "lib_alt_Accel/*.cpp"
        "lib_alt_Accel/Eigen/*.cpp"
        "lib_alt_Accel/Field/*.cpp"
        "lib_alt_Accel/Fopr/*.cpp"
        "lib_alt_Accel/Measurements/Fermion/*.cpp"
        "lib_alt_Accel/Measurements/Gauge/*.cpp"
        "lib_alt_Accel/Solver/*.cpp"
    )
    # Only compile the deviceManager that matches the selected backend
    if(USE_ACCEL_IMPL STREQUAL "cuda")
        list(APPEND BRIDGE_ACCEL_SOURCES
            "${CMAKE_CURRENT_SOURCE_DIR}/lib_alt_Accel/ResourceManager/deviceManager_CUDA.cpp")
    elseif(USE_ACCEL_IMPL STREQUAL "nvshmem")
        list(APPEND BRIDGE_ACCEL_SOURCES
            "${CMAKE_CURRENT_SOURCE_DIR}/lib_alt_Accel/ResourceManager/deviceManager_NVSHMEM.cpp")
    elseif(USE_ACCEL_IMPL STREQUAL "openacc")
        list(APPEND BRIDGE_ACCEL_SOURCES
            "${CMAKE_CURRENT_SOURCE_DIR}/lib_alt_Accel/ResourceManager/deviceManager_OpenACC.cpp")
        file(GLOB ACC_BRIDGE CONFIGURE_DEPENDS "lib_alt_Accel/BridgeOpenACC/*.cpp")
        list(APPEND BRIDGE_ACCEL_SOURCES ${ACC_BRIDGE})
    endif()
endif()

# lib_alt_QXS: Fujitsu A64FX SIMD (Fugaku)
set(BRIDGE_QXS_SOURCES)
if(USE_ALT_QXS)
    file(GLOB BRIDGE_QXS_SOURCES CONFIGURE_DEPENDS
        "lib_alt_QXS/*.cpp"
        "lib_alt_QXS/Action/Fermion/*.cpp"
        "lib_alt_QXS/BridgeQXS/*.cpp"
        "lib_alt_QXS/Eigen/*.cpp"
        "lib_alt_QXS/extra/*.cpp"
        "lib_alt_QXS/Field/*.cpp"
        "lib_alt_QXS/Fopr/*.cpp"
        "lib_alt_QXS/Force/Fermion/*.cpp"
        "lib_alt_QXS/Measurements/Fermion/*.cpp"
        "lib_alt_QXS/Measurements/Gauge/*.cpp"
        "lib_alt_QXS/Smear/*.cpp"
        "lib_alt_QXS/Solver/*.cpp"
    )
endif()

# ----------- target -----------

add_library(bridge STATIC
    ${BRIDGE_CORE}
    ${BRIDGE_COMMUN}
    ${BRIDGE_OPT}
    ${BRIDGE_THREAD}
    ${BRIDGE_ALT_SOURCES}
    ${BRIDGE_ACCEL_SOURCES}
    ${BRIDGE_QXS_SOURCES}
)

# ----------- include paths (PUBLIC → propagated to every consumer) -----------

target_include_directories(bridge PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/lib"
)
if(USE_ALTERNATIVE)
    target_include_directories(bridge PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/lib_alt")
endif()
if(USE_ALT_ACCEL)
    target_include_directories(bridge PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/lib_alt_Accel")
endif()
if(USE_ALT_QXS)
    target_include_directories(bridge PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/lib_alt_QXS")
endif()

# ----------- compile definitions (PUBLIC → propagated) -----------

target_compile_definitions(bridge PUBLIC
    USE_FACTORY
    USE_STD_COMPLEX
    LIB_CPP11
    USE_GROUP_${BRIDGE_GAUGE_GROUP}
    $<$<CONFIG:Debug>:DEBUG>
    $<$<NOT:$<CONFIG:Debug>>:NDEBUG>
)

if(USE_OPT_CODE)
    target_compile_definitions(bridge PUBLIC USE_IMP)
else()
    target_compile_definitions(bridge PUBLIC USE_ORG)
endif()

if(USE_MPI STREQUAL "fjmpi")
    target_compile_definitions(bridge PUBLIC USE_MPI USE_FJMPI)
elseif(USE_MPI)
    target_compile_definitions(bridge PUBLIC USE_MPI)
endif()

if(USE_OPENMP)
    target_compile_definitions(bridge PUBLIC USE_OPENMP)
endif()

if(USE_ALTERNATIVE)
    target_compile_definitions(bridge PUBLIC USE_ALT_CODE)
endif()

if(USE_ALT_QXS)
    target_compile_definitions(bridge PUBLIC USE_ALT_QXS)
    if(USE_QXS_ARCH STREQUAL "acle")
        target_compile_definitions(bridge PUBLIC USE_QXS_ACLE)
    endif()
endif()

if(USE_ALT_ACCEL)
    target_compile_definitions(bridge PUBLIC USE_ALT_ACCEL)
    if(USE_ACCEL_IMPL STREQUAL "cuda")
        target_compile_definitions(bridge PUBLIC USE_ACCEL_CUDA)
        if(USE_DIRECT_COMM)
            target_compile_definitions(bridge PUBLIC USE_DIRECT_COMM)
        endif()
    elseif(USE_ACCEL_IMPL STREQUAL "nvshmem")
        target_compile_definitions(bridge PUBLIC USE_ACCEL_NVSHMEM)
    elseif(USE_ACCEL_IMPL STREQUAL "openacc")
        target_compile_definitions(bridge PUBLIC USE_ACCEL_OPENACC)
    endif()
endif()

# nvcc/nvc++ don't provide strdup in <string.h> by default
if(CMAKE_CXX_COMPILER_ID MATCHES "NVHPC|PGI")
    target_compile_definitions(bridge PUBLIC USE_STRDUP)
endif()

# ----------- external dependencies -----------

if(USE_MPI STREQUAL "fjmpi")
    # FJMPI is integrated into the mpiFCCpx compiler wrapper; no find_package needed
elseif(USE_MPI)
    find_package(MPI REQUIRED CXX)
    target_link_libraries(bridge PUBLIC MPI::MPI_CXX)
endif()

if(USE_OPENMP)
    find_package(OpenMP REQUIRED CXX)
    target_link_libraries(bridge PUBLIC OpenMP::OpenMP_CXX)
endif()

if(USE_ALT_ACCEL AND (USE_ACCEL_IMPL STREQUAL "cuda" OR USE_ACCEL_IMPL STREQUAL "nvshmem"))
    # CUDAToolkit provides cuda_runtime.h path needed by channel.h and deviceManager_*.cpp
    find_package(CUDAToolkit REQUIRED)
    target_include_directories(bridge PUBLIC ${CUDAToolkit_INCLUDE_DIRS})
    target_link_libraries(bridge PUBLIC bridgecuda CUDA::cudart)
    if(USE_ACCEL_IMPL STREQUAL "nvshmem")
        find_library(NVSHMEM_LIB nvshmem
            HINTS ENV NVSHMEM_HOME
            PATH_SUFFIXES lib lib64)
        if(NVSHMEM_LIB)
            target_link_libraries(bridge PUBLIC ${NVSHMEM_LIB})
        endif()
    endif()
endif()

if(USE_ALT_ACCEL AND USE_ACCEL_IMPL STREQUAL "openacc")
    find_package(OpenACC REQUIRED CXX)
    # OpenACC::OpenACC_CXX propagates -acc (nvc++) or -fopenacc (GCC) as both compile and link flags
    target_link_libraries(bridge PUBLIC OpenACC::OpenACC_CXX)
endif()

# ----------- optional external libraries -----------

if(USE_LIME)
    find_path(LIME_INCLUDE_DIR limeWriter.h
        HINTS "${CMAKE_SOURCE_DIR}/extra/lime-1.3.2/include")
    find_library(LIME_LIBRARY lime
        HINTS "${CMAKE_SOURCE_DIR}/extra/lime-1.3.2/lib")
    target_include_directories(bridge PUBLIC ${LIME_INCLUDE_DIR})
    target_link_libraries(bridge PUBLIC ${LIME_LIBRARY})
    target_compile_definitions(bridge PUBLIC USE_LIMELIB)
endif()

if(USE_SFMT)
    find_path(SFMT_INCLUDE_DIR SFMT.h
        HINTS "${CMAKE_SOURCE_DIR}/extra/SFMT-1.5.1_bridge")
    find_library(SFMT_LIBRARY sfmt
        HINTS "${CMAKE_SOURCE_DIR}/extra/SFMT-1.5.1_bridge")
    target_include_directories(bridge PUBLIC ${SFMT_INCLUDE_DIR})
    target_link_libraries(bridge PUBLIC ${SFMT_LIBRARY})
    target_compile_definitions(bridge PUBLIC USE_SFMTLIB)
endif()

if(USE_NTL)
    find_path(NTL_INCLUDE_DIR NTL/ZZ.h
        HINTS "${CMAKE_SOURCE_DIR}/extra/ntl/include")
    find_library(NTL_LIBRARY ntl
        HINTS "${CMAKE_SOURCE_DIR}/extra/ntl/lib")
    find_library(GMP_LIBRARY gmp)
    target_include_directories(bridge PUBLIC ${NTL_INCLUDE_DIR})
    target_link_libraries(bridge PUBLIC ${NTL_LIBRARY} ${GMP_LIBRARY})
    target_compile_definitions(bridge PUBLIC USE_NTLLIB)
endif()

if(USE_FFTW)
    find_path(FFTW_INCLUDE_DIR fftw3.h
        HINTS "${CMAKE_SOURCE_DIR}/extra/fftw/include")
    set(_FFTW_LIBS)
    if(USE_MPI)
        find_library(FFTW_MPI_LIBRARY fftw3_mpi
            HINTS "${CMAKE_SOURCE_DIR}/extra/fftw/lib")
        list(APPEND _FFTW_LIBS ${FFTW_MPI_LIBRARY})
    endif()
    if(USE_OPENMP)
        find_library(FFTW_OMP_LIBRARY fftw3_omp
            HINTS "${CMAKE_SOURCE_DIR}/extra/fftw/lib")
        list(APPEND _FFTW_LIBS ${FFTW_OMP_LIBRARY})
    endif()
    find_library(FFTW_LIBRARY fftw3
        HINTS "${CMAKE_SOURCE_DIR}/extra/fftw/lib")
    list(APPEND _FFTW_LIBS ${FFTW_LIBRARY})
    target_include_directories(bridge PUBLIC ${FFTW_INCLUDE_DIR})
    target_link_libraries(bridge PUBLIC ${_FFTW_LIBS})
    target_compile_definitions(bridge PUBLIC USE_FFTWLIB)
endif()

if(USE_YAMLCPP)
    find_path(YAMLCPP_INCLUDE_DIR yaml-cpp/yaml.h
        HINTS "${CMAKE_SOURCE_DIR}/extra/yaml-cpp/include")
    find_library(YAMLCPP_LIBRARY yaml-cpp
        HINTS "${CMAKE_SOURCE_DIR}/extra/yaml-cpp/lib")
    target_include_directories(bridge PUBLIC ${YAMLCPP_INCLUDE_DIR})
    target_link_libraries(bridge PUBLIC ${YAMLCPP_LIBRARY})
    target_compile_definitions(bridge PUBLIC USE_YAMLCPPLIB)
endif()

if(USE_TINYXML2)
    find_path(TINYXML2_INCLUDE_DIR tinyxml2.h
        HINTS "${CMAKE_SOURCE_DIR}/extra/tinyxml2/include")
    find_library(TINYXML2_LIBRARY tinyxml2
        HINTS "${CMAKE_SOURCE_DIR}/extra/tinyxml2/lib")
    target_include_directories(bridge PUBLIC ${TINYXML2_INCLUDE_DIR})
    target_link_libraries(bridge PUBLIC ${TINYXML2_LIBRARY})
    target_compile_definitions(bridge PUBLIC USE_XML USE_TINYXML2LIB)
endif()

if(USE_QWS)
    find_path(QWS_INCLUDE_DIR qws.h
        HINTS "${CMAKE_SOURCE_DIR}/extra/qws-modified")
    find_library(QWS_LIBRARY qws
        HINTS "${CMAKE_SOURCE_DIR}/extra/qws-modified")
    target_include_directories(bridge PUBLIC ${QWS_INCLUDE_DIR})
    target_link_libraries(bridge PUBLIC ${QWS_LIBRARY})
    target_compile_definitions(bridge PUBLIC USE_QWSLIB)
endif()
