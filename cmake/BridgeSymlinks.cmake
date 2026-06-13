# BridgeSymlinks.cmake
# ---------------------------------------------------------------------------
# Recreate the source-tree symlinks the Makefile build's `link-inline` target
# made, so #include "lib_alt_Accel/BridgeACC/..." and "lib_alt_QXS/inline/..."
# resolve to the backend/arch selected by the cache.  Included from the root
# CMakeLists.txt (runs in its scope).  No targets, only filesystem setup.

# ----------- BridgeACC symlink -----------
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

# ----------- QXS inline symlink -----------
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
