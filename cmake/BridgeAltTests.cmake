# BridgeAltTests.cmake
# ---------------------------------------------------------------------------
# Build the standalone test_alt_* programs (flat *.cpp-in-root layout, or the
# src/ + src/Tests/ layout) and link them against the bridge library.  Included
# from the root CMakeLists.txt AFTER add_subdirectory(src)/add_subdirectory(tests)
# so `bridge` (and bridgecuda_mrhs, for the batched-AMG test) already exist.

file(READ "${CMAKE_SOURCE_DIR}/VERSION" _BRIDGE_VERSION)
string(STRIP "${_BRIDGE_VERSION}" _BRIDGE_VERSION)

foreach(_dir IN ITEMS
        test_alt_multigrid
        test_alt_dw_multigrid
        test_alt_dw_multigrid_dev2
        test_alt_dw_multigrid_dev4
        test_alt_dw_multigrid_dev5
        test_alt_dw_multigrid_dev6
        test_alt_amg_dw_spmv
        test_alt_amg_dw_GEMM
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

# ----------- amg_dw_GEMM (was dev8): link the batched AMG MRHS kernels ---------
# The MRHS / tensor-core kernels live in the additive bridgecuda_mrhs static lib
# (src/lib_alt_Accel/BridgeCUDA/bridgeACC_mrhs_block_float.cu) instead of being
# compiled into this exe.  Linking it pulls in the mrhs_live:: API and its PUBLIC
# include dir (mrhs_block_live.h), so the test stays a plain C++ target.
if(TARGET bridge_test_alt_amg_dw_GEMM AND TARGET bridgecuda_mrhs)
    target_link_libraries(bridge_test_alt_amg_dw_GEMM PRIVATE bridgecuda_mrhs)
endif()
