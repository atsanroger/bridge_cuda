################################################################################
#        @author  Hideo Matsufuru (matufuru) 
#                 $LastChangedBy: wlchen $
#        @date    $LastChangedDate:: 2024-10-02 19:34:22 #$
#        @version $LastChangedRevision: 2562 $
################################################################################
SHELL := /bin/bash

ARCH = -gencode=arch=compute_80,code=sm_80 \
       -gencode=arch=compute_86,code=sm_86 \
       -gencode=arch=compute_90,code=sm_90

sources = bridgeACC.cu \
          bridgeACC_double.cu \
          bridgeACC_float.cu \
          bridgeACC_AField.cu \
          bridgeACC_AField_double.cu \
          bridgeACC_AField_float.cu \
          bridgeACC_Wilson_double.cu \
          bridgeACC_Wilson_float.cu \
          bridgeACC_Clover_double.cu \
          bridgeACC_Clover_float.cu \
          bridgeACC_CloverTerm_double.cu \
          bridgeACC_CloverTerm_float.cu \
          bridgeACC_Clover_coarse_double.cu \
          bridgeACC_Clover_coarse_float.cu \
          constant_memory.cu \
          bridgeACC_Domainwall_5din_double.cu \
          bridgeACC_Domainwall_5din_float.cu

objects = $(sources:%.cu=%.o)

CXX = nvcc
CXXFLAGS = -O3 -std=c++11 $(ARCH) -rdc=true --disable-warnings
CXXFLAGS += -DUSE_GROUP_SU3 $(CFLAGS)

LDLIBS = -lcudart

RANLIB ?= ranlib

TOTAL := $(words $(sources))
COUNT_FILE := .progress_count

define UPDATE_PROGRESS
	@COUNT=$$(($(shell cat $(COUNT_FILE)) + 1)); \
	echo $$COUNT > $(COUNT_FILE); \
	printf "\rProgress: [$$COUNT/$(TOTAL)] Compiling $1, "
endef

all : start_message libbridgecuda.a

start_message:
	@echo "Bridge++'s CUDA LIBRARY STARTING"
	@echo "Total source files: $(TOTAL)"

libbridgecuda.a : $(objects)
	$(CXX) -dlink $(objects) -o link.o $(ARCH) $(LDLIBS)
	ar rcs $@ $(objects) link.o
	$(RANLIB) $@
	@echo "libbridgecuda.a has been successfully compiled."

$(objects): $(COUNT_FILE)

$(COUNT_FILE):
	@echo 0 > $(COUNT_FILE)

%.o : %.cu
	$(call UPDATE_PROGRESS,$<)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f libbridgecuda.a $(objects) link.o *.a *.so *.elf $(COUNT_FILE) *.o
	find . -name '*.o' -delete
	find . -name '*.a' -delete
	find . -name '*.so' -delete
	@echo "libbridgecuda.a and all object files have been successfully cleaned."