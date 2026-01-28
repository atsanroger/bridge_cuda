ARCH = \
 -gencode=arch=compute_80,code=sm_80 \
 -gencode=arch=compute_90,code=sm_90

# -gencode=arch=compute_86,code=sm_86 \
#-gencode=arch=compute_61,code=sm_61 \
# -gencode=arch=compute_70,code=sm_70 \
# -gencode=arch=compute_100,code=sm_100

sources = bridgeACC.cu \
 bridgeACC_double.cu \
 bridgeACC_float.cu \
 bridgeACC_AField.cu \
 bridgeACC_AField_double.cu \
 bridgeACC_AField_float.cu  \
 bridgeACC_Wilson_double.cu \
 bridgeACC_Wilson_float.cu  \
 bridgeACC_Clover_double.cu \
 bridgeACC_Clover_float.cu  \
 bridgeACC_CloverTerm_double.cu \
 bridgeACC_CloverTerm_float.cu  \
 bridgeACC_Clover_coarse_double.cu \
 bridgeACC_Clover_coarse_float.cu \
 bridgeACC_Domainwall_5din_double.cu \
 bridgeACC_Domainwall_5din_float.cu \
# bridgeACC_Domainwall_coarse_double.cu \
# bridgeACC_Domainwall_coarse_float.cu \

objects = $(sources:%.cu=%.o)

NVSHMEM_HOME = /work/opt/local/x86_64/cores/nvidia_gcc/23.3_12.2.0/Linux_x86_64/23.3/comm_libs/nvshmem


CXX = nvcc
CXXFLAGS = -O3 $(ARCH) -rdc=true --std=c++11
CXXFLAGS += -DUSE_GROUP_SU3 -I$(NVSHMEM_HOME)/include

LDLIBS = -L$(NVSHMEM_HOME)/lib -lm -lnvshmem -lcudadevrt -lcudart

AR = ar r
RANLIB = ranlib

all : libbridgecuda.a

libbridgecuda.a : $(objects)
	$(CXX) -dlink $(objects) -o link.o $(LDLIBS) $(ARCH)
	$(AR) $@ $(objects) link.o
	$(RANLIB) $@

%.o : %.cu
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f libbridgecuda.a $(objects) link.o
