# ============================================================
#  Bridge++ test suite  →  bridge_test
# ============================================================

option(USE_TESTMANAGER "Enable test manager (auto-register)" OFF)

# ----------- test module list -----------
set(_TEST_MODULES
    Eigensolver
    FFT
    Gauge
    GradientFlow
    HMC
    HotStart
    IO
    Mult
    PolyakovLoop
    QuarkNumberSusceptibility
    RandomNumbers
    Rational
    SF_fAfP
    Solver
    Spectrum
    TopologicalCharge
    WilsonLoop
)

# ----------- source collection -----------
file(GLOB _TEST_SRCS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

foreach(_mod IN LISTS _TEST_MODULES)
    file(GLOB _mod_srcs CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/src/${_mod}/*.cpp")
    list(APPEND _TEST_SRCS ${_mod_srcs})
endforeach()

# ----------- executable -----------
add_executable(bridge_test ${_TEST_SRCS})

target_include_directories(bridge_test PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
)

target_link_libraries(bridge_test PRIVATE bridge)

# ----------- compile definitions (USE_TEST_<MODULE>) -----------
# QuarkNumberSusceptibility  →  QUARKNUMSUSCEPT  (Makefile substitution)
# SF_fAfP  →  SF_FAFP  (tr a-z A-Z)
foreach(_mod IN LISTS _TEST_MODULES)
    string(TOUPPER "${_mod}" _mod_upper)
    string(REPLACE "QUARKNUMBERSUSCEPTIBILITY" "QUARKNUMSUSCEPT" _mod_upper "${_mod_upper}")
    target_compile_definitions(bridge_test PRIVATE "USE_TEST_${_mod_upper}")
endforeach()

if(USE_TESTMANAGER)
    target_compile_definitions(bridge_test PRIVATE
        USE_TESTMANAGER
        USE_TESTMANAGER_AUTOREGISTER
    )
endif()

# ----------- copy ParameterFiles to binary dir at configure time -----------

set(_PARAM_DIR "${CMAKE_CURRENT_SOURCE_DIR}/ParameterFiles")
set(_BIN_DIR   "${CMAKE_CURRENT_BINARY_DIR}")

# main.yaml
configure_file("${_PARAM_DIR}/main.yaml" "${_BIN_DIR}/main.yaml" COPYONLY)

# ParameterFiles/<Module>/test_*.yaml  (one glob per module)
foreach(_mod IN LISTS _TEST_MODULES)
    file(GLOB _yaml_files "${_PARAM_DIR}/${_mod}/test_*.yaml")
    foreach(_f IN LISTS _yaml_files)
        get_filename_component(_fname "${_f}" NAME)
        configure_file("${_f}" "${_BIN_DIR}/${_fname}" COPYONLY)
    endforeach()
endforeach()

# parameter_rational* files
file(GLOB _rational_files "${_PARAM_DIR}/parameter_rational*")
foreach(_f IN LISTS _rational_files)
    get_filename_component(_fname "${_f}" NAME)
    configure_file("${_f}" "${_BIN_DIR}/${_fname}" COPYONLY)
endforeach()

# ----------- decompress sample configuration at configure time -----------

set(_CONF_GZ "${CMAKE_CURRENT_SOURCE_DIR}/Configurations/conf_04040408.txt.gz")
set(_CONF_OUT "${_BIN_DIR}/conf_04040408.txt")

if(EXISTS "${_CONF_GZ}" AND NOT EXISTS "${_CONF_OUT}")
    execute_process(
        COMMAND gunzip -c "${_CONF_GZ}"
        OUTPUT_FILE "${_CONF_OUT}"
    )
endif()
