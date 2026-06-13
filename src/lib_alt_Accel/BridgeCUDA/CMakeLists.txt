set(BRIDGECUDA_SOURCES
    bridgeACC.cu
    bridgeACC_double.cu
    bridgeACC_float.cu
    bridgeACC_AField.cu
    bridgeACC_AField_double.cu
    bridgeACC_AField_float.cu
    bridgeACC_Wilson_double.cu
    bridgeACC_Wilson_float.cu
    bridgeACC_Clover_double.cu
    bridgeACC_Clover_float.cu
    bridgeACC_CloverTerm_double.cu
    bridgeACC_CloverTerm_float.cu
    bridgeACC_Clover_coarse_double.cu
    bridgeACC_Clover_coarse_float.cu
    constant_memory.cu
    bridgeACC_Domainwall_5din_double.cu
    bridgeACC_Domainwall_5din_float.cu
    bridgeACC_Domainwall_coarse_double.cu
    bridgeACC_Domainwall_coarse_float.cu
)

add_library(bridgecuda STATIC ${BRIDGECUDA_SOURCES})

# USE_GROUP_SU3 / USE_GROUP_SU2 / USE_GROUP_SU_N  (from BRIDGE_GAUGE_GROUP)
target_compile_definitions(bridgecuda PRIVATE USE_GROUP_${BRIDGE_GAUGE_GROUP})

target_compile_options(bridgecuda PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:--disable-warnings>
)

set_target_properties(bridgecuda PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON  # device-link embedded into the archive, matching Makefile behaviour
    CUDA_STANDARD              11
    CUDA_STANDARD_REQUIRED     ON
)

# ---------------------------------------------------------------------------
# bridgecuda_mrhs: AMG-on-A batched MRHS / tensor-core kernels (TF32 mma.sync).
# ADDITIVE separate static lib so the existing bridgecuda (C++11) and every
# shared fermion/Corr2pt/CGNE path stay byte-for-byte untouched.  This TU needs
# C++17 + USE_TENSORCORE_MRHS, which is why it is its own target rather than
# folded into bridgecuda.  It exposes the mrhs_live:: API (mrhs_block_live.h);
# any test/solver that wants the batched AMG kernels just links bridgecuda_mrhs.
# (Was: compiled into the dev8 exe.)  CUDA_ARCHITECTURES is intentionally NOT set
# -- it inherits the global CMAKE_CUDA_ARCHITECTURES from the active preset
# (cuda-ampere-mpi=80;86, cuda-hopper-mpi=90, ...), exactly like bridgecuda; the
# m16n8k8 TF32 mma.sync is available on every sm_80+ target the presets build.
add_library(bridgecuda_mrhs STATIC
    ${CMAKE_CURRENT_SOURCE_DIR}/bridgeACC_mrhs_block_float.cu)

target_include_directories(bridgecuda_mrhs PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR})                       # bridgeACC.h, mrhs_block_live.h, src/*-inc.h
target_include_directories(bridgecuda_mrhs PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/lib
    ${CMAKE_SOURCE_DIR}/src/lib_alt
    ${CMAKE_SOURCE_DIR}/src/lib_alt_Accel)

target_compile_definitions(bridgecuda_mrhs PRIVATE
    USE_GROUP_${BRIDGE_GAUGE_GROUP}
    USE_TENSORCORE_MRHS)

target_compile_options(bridgecuda_mrhs PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:--disable-warnings>)

set_target_properties(bridgecuda_mrhs PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON
    CUDA_STANDARD              17
    CUDA_STANDARD_REQUIRED     ON)

# host-side symbols (dev_ptr, ...) live in bridgecuda; device code here is
# self-contained (module-local Moebius b/c), so this is a host-link dependency.
target_link_libraries(bridgecuda_mrhs PUBLIC bridgecuda)
