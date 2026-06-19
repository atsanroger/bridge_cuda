.SUFFIXES: .cpp .h .o

BRIDGE_INSTALL_PATH = ../../build

include $(BRIDGE_INSTALL_PATH)/make.inc

headers =
sources = pyb.cpp \
	pyb_common.cpp \
	pyb_params.cpp \
	pyb_tools_su_n.cpp \
	pyb_field.cpp \
	pyb_index.cpp \
	pyb_fopr.cpp \
	pyb_solver.cpp \
	pyb_smear.cpp \
	pyb_shift_field.cpp \
	pyb_gauge_fixing.cpp \
	pyb_wilsonloop.cpp \
	pyb_source.cpp \
	pyb_fprop.cpp \
	pyb_gammamatrix.cpp \
	pyb_correlator.cpp \
	pyb_staple.cpp \
	pyb_action.cpp \
	pyb_force.cpp \
	pyb_hmc.cpp \
	pyb_rng.cpp \
	pyb_fft.cpp \
	pyb_gradientflow.cpp \
	pyb_topo.cpp \
	pyb_qsuscept.cpp \
	pyb_eigen.cpp \
	pyb_io.cpp \
	pyb_com.cpp \
	pyb_tools.cpp


objects = $(sources:%.cpp=%.o)

program = pybridge$(shell python3-config --extension-suffix)

# python paths
#
PYTHON_INCLUDES := $(shell python3-config --includes)

# pybind11 paths
#
PYBIND11_INCLUDES := -I../pybind11/include

# mpi4py paths
#
ifeq ($(findstring -DUSE_MPI,$(CPPFLAGS)),-DUSE_MPI)
  MPI4PY_INCLUDES := $(shell python3 -c "import mpi4py; print('-DUSE_MPI4PY -I'+mpi4py.get_include())")
endif


CXXFLAGS += $(PYTHON_INCLUDES) $(PYBIND11_INCLUDES) $(MPI4PY_INCLUDES)
CPPFLAGS += -DTRUNK

all: msg $(program)

$(program): $(objects)
	$(LD) -shared $(LDFLAGS) -o $@ $(objects) $(LDLIBS)

%.o: %.cpp $(headers)
	$(CXX) -fPIC $(CXXFLAGS) -c -o $@ $<

msg:
	@cat $(BRIDGE_INSTALL_PATH)/message.txt

clean:
	rm -f $(program) $(objects)

