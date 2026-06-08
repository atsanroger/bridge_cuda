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
