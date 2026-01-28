#####  GNU C/C++  ##############################################

ifeq ($(strip $(target)),PC_GNU)

ifeq ($(use_mpi),yes)
  CXX = mpic++
else
  CXX = g++
endif

ifeq ($(use_thread),omp)
  CXXFLAGS += -fopenmp
  LDFLAGS  += -fopenmp
endif

ifeq ($(use_opt_level),low)
  OPTFLAGS = -O0
endif
ifeq ($(use_opt_level),std)
  OPTFLAGS = -O
endif
ifeq ($(use_opt_level),high)
  OPTFLAGS = -O3
endif

ifeq ($(use_static),yes)
  LDFLAGS  += -static
endif

ifeq ($(use_debug),yes)
  CXXFLAGS += -g -Wall -Wno-unused -Wno-sign-compare -Wno-unknown-pragmas -fbounds-check -fsanitize=undefined -fno-sanitize-recover
  CPP_DEFS += -DDEBUG
  LDFLAGS += -fsanitize=undefined -fno-sanitize-recover
else
  CPP_DEFS += -DNDEBUG
endif

ifeq ($(use_cpp11),yes)
  CXXFLAGS += -std=gnu++11
endif

#CXXFLAGS += -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64
CXXFLAGS += $(OPTFLAGS)

ifeq ($(use_alt_simd),yes)
  ifeq ($(simd_arch),AVX1)
    CXXFLAGS += -mavx
  endif
  ifeq ($(simd_arch),AVX2)
    CXXFLAGS += -mavx2 -mfma
  endif
  ifeq ($(simd_arch),AVX512)
    CXXFLAGS += -mavx512f
  endif
endif

LD = $(CXX)
LDLIBS = -lm

CXXDEP = $(CXX)
CXXDEPFLAGS = -MM

AR = ar r
RANLIB = ranlib

TARGET_MACHINE = PC_GNU

#####  CLANG C/C++  ##############################################

else ifeq ($(strip $(target)),PC_CLANG)

ifeq ($(use_mpi),yes)
  CXX = mpic++
else
  CXX = clang++
endif

ifeq ($(use_thread),omp)
  CXXFLAGS += -fopenmp
  LDFLAGS  += -fopenmp
endif

ifeq ($(use_opt_level),low)
  OPTFLAGS = -O0
endif
ifeq ($(use_opt_level),std)
  OPTFLAGS = -O
endif
ifeq ($(use_opt_level),high)
  OPTFLAGS = -O3
endif

ifeq ($(use_static),yes)
  LDFLAGS  += -static
endif

ifeq ($(use_debug),yes)
  CXXFLAGS += -g -Weverything
  CPP_DEFS += -DDEBUG
  LDFLAGS += -fsanitize=undefined -fno-sanitize-recover
else
  CPP_DEFS += -DNDEBUG
endif

ifeq ($(use_cpp11),yes)
  CXXFLAGS += -std=gnu++11
endif

ifeq ($(use_python),yes)
  CXXFLAGS += -fPIC
endif

#CXXFLAGS += -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64
CXXFLAGS += $(OPTFLAGS)

LD = $(CXX)
LDLIBS = -lm

CXXDEP = $(CXX)
CXXDEPFLAGS = -MM

AR = ar r
RANLIB = ranlib

TARGET_MACHINE = PC_CLANG

#####  Intel C++  ##############################################

else ifeq ($(strip $(target)),PC_INTEL)

ifeq ($(use_mpi),yes)
  #CXX = mpic++
  CXX = mpiicpc
  CXXFLAGS += -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX
else
  CXX = icpc
endif

ifeq ($(use_thread),omp)
  CXXFLAGS += -qopenmp  # -openmp
  LDFLAGS  += -qopenmp  # -openmp
endif

ifeq ($(use_opt_level),low)
  OPTFLAGS = -O0
endif
ifeq ($(use_opt_level),std)
  OPTFLAGS = -O2
endif
ifeq ($(use_opt_level),high)
  #OPTFLAGS = -fast
  OPTFLAGS = -O3 -no-prec-div -xHost  # -ipo
endif

ifeq ($(use_static),yes)
  LDFLAGS  += -static-intel  # -static
  ## for intel MPI library
  # LDFLAGS  += -static_mpi
else
  CXXFLAGS += -mcmodel=medium  # -shared-intel
  LDFLAGS  += -mcmodel=medium  # -shared-intel
endif

ifeq ($(use_debug),yes)
  CXXFLAGS += -g -debug all -traceback
  CPP_DEFS += -DDEBUG
else
  CPP_DEFS += -DNDEBUG
endif

ifeq ($(use_cpp11),yes)
  CXXFLAGS += -std=c++11
  ## NB. cpp11+c99 causes compile errors for icpc 15.0.3
  #  $(error Incompatible options: cpp11 conflicts c99 for INTEL)
  $(warning cpp11 may fail for PC_INTEL, try cpp11=no)
endif

CXXFLAGS += $(OPTFLAGS)

LD = $(CXX)
LDLIBS = -lm

CXXDEP = $(CXX)
CXXDEPFLAGS = -MM

AR = ar r
RANLIB = ranlib

TARGET_MACHINE = PC_INTEL

#####  NVIDIA HPC SDK C++  #####################################

else ifeq ($(strip $(target)),PC_NVIDIA)

CPP_DEFS += -DUSE_STRDUP

ifeq ($(use_mpi),yes)
  CXX = mpic++
else
  CXX = nvc++
endif

ifeq ($(use_thread),omp)
  CXXFLAGS += -mp
  LDFLAGS  += -mp
endif

ifeq ($(use_opt_level),low)
  OPTFLAGS = -O0
endif
ifeq ($(use_opt_level),std)
  OPTFLAGS = -O2
endif
ifeq ($(use_opt_level),high)
  OPTFLAGS = -fast -O3
endif

ifeq ($(use_static),yes)
  LDFLAGS  += -Wl,-Bstatic
endif

ifeq ($(use_debug),yes)
  CXXFLAGS += -g -Minform=inform -Mbounds -Ktrap=fp
  CPP_DEFS += -DDEBUG
else
  CPP_DEFS += -DNDEBUG
endif

ifeq ($(use_cpp11),yes)
  CXXFLAGS += -std=c++11
endif

CXXFLAGS += $(OPTFLAGS)

LD = $(CXX)
LDLIBS = -lm

ifeq ($(use_fftw_library),yes)
  ifeq ($(use_thread),omp)
    LDLIBS += -lgomp
  endif
endif

ifeq ($(use_alt_accel),yes)

ifeq ($(use_accel_impl),openacc)
  GPU_SUFFIX :=
#  GPU_SUFFIX := mem:separate:pinnedalloc
#  GPU_SUFFIX := pinned
  CXXFLAGS += -acc
  LDFLAGS  += -acc
else ifeq ($(use_accel_impl),cuda)
  GPU_SUFFIX :=
endif

ifeq ($(use_device_arch),kepler)
  CXXFLAGS += -gpu=ptxinfo,cc35,$(GPU_SUFFIX),maxregcount:
  LDFLAGS  += -gpu=ptxinfo,cc35,$(GPU_SUFFIX),maxregcount:
endif
ifeq ($(use_device_arch),pascal)
  CXXFLAGS += -gpu=ptxinfo,cc60,$(GPU_SUFFIX),maxregcount:
  LDFLAGS  += -gpu=ptxinfo,cc60,$(GPU_SUFFIX),maxregcount:
endif
ifeq ($(use_device_arch),volta)
  CXXFLAGS += -gpu=ptxinfo,cc70,$(GPU_SUFFIX),maxregcount:
  LDFLAGS  += -gpu=ptxinfo,cc70,$(GPU_SUFFIX),maxregcount:
endif
ifeq ($(use_device_arch),ampere)
  CXXFLAGS += -gpu=ptxinfo,cc80,$(GPU_SUFFIX),maxregcount:
  LDFLAGS  += -gpu=ptxinfo,cc80,$(GPU_SUFFIX),maxregcount:
endif
ifeq ($(use_device_arch),hopper)
  CXXFLAGS += -gpu=ptxinfo,cc90,$(GPU_SUFFIX),maxregcount:
  LDFLAGS  += -gpu=ptxinfo,cc90,$(GPU_SUFFIX),maxregcount:
endif

ifeq ($(use_device_arch),rumors)
  CXXFLAGS += -gpu=ptxinfo,cc100,$(GPU_SUFFIX),maxregcount:
  LDFLAGS  += -gpu=ptxinfo,cc100,$(GPU_SUFFIX),maxregcount:
endif

ifeq ($(use_device_arch),auto)
  # Use nvidia-smi to detect the GPU architecture
  GPU_ARCH := $(shell nvidia-smi --query-gpu=compute_cap --format=csv,noheader,nounits | awk -F. '{print "cc" $$1 $$2}')
  CXXFLAGS += -gpu=ptxinfo,$(GPU_ARCH),maxregcount:
  LDFLAGS  += -gpu=ptxinfo,$(GPU_ARCH),maxregcount:
endif

ifeq ($(use_accel_impl),cuda)
  CXXFLAGS += -cuda
  LDFLAGS  += -cuda
endif

endif # use_alt_accel


CXXDEP = $(CXX)
CXXDEPFLAGS = -MM

AR = ar r
RANLIB = ranlib

TARGET_MACHINE = PC_NVIDIA


#####  NVIDIA CUDA  ############################################

else ifeq ($(strip $(target)),CUDA)

CPP_DEFS += -DUSE_STRDUP

ifeq ($(use_mpi),yes)
  CXX = nvcc -ccbin mpic++ 
else
  CXX = nvcc
endif

ifeq ($(use_thread),omp)
  CXXFLAGS += -Xcompiler -fopenmp
  LDFLAGS  += -Xcompiler -fopenmp
endif

ifeq ($(use_opt_level),low)
  OPTFLAGS = -O0
endif
ifeq ($(use_opt_level),std)
  OPTFLAGS = -O2
endif
ifeq ($(use_opt_level),high)
  OPTFLAGS = -O3
endif

ifeq ($(use_static),yes)
  LDFLAGS  += -ccbin mpic++
endif

ifeq ($(use_debug),yes)
  CXXFLAGS += -g
  CPP_DEFS += -DDEBUG
else
  CPP_DEFS += -DNDEBUG
endif

ifeq ($(use_cpp11),yes)
  CXXFLAGS += -std=c++11
endif

CXXFLAGS += $(OPTFLAGS)

LD = $(CXX)
LDLIBS = -lm -lnvidia-ml -lcuda -lcudart

ifeq ($(use_fftw_library),yes)
  ifeq ($(use_thread),omp)
    LDLIBS += -lgomp
  endif
endif

ifeq ($(use_alt_accel),yes)

ifeq ($(use_device_arch),kepler)
endif
ifeq ($(use_device_arch),pascal)
  CXXFLAGS += -gpu=ptxinfo,cc60,maxregcount:
  LDFLAGS  += -gpu=ptxinfo,cc60,maxregcount:
endif
ifeq ($(use_device_arch),volta)
  CXXFLAGS += -gpu=ptxinfo,cc70,maxregcount:
  LDFLAGS  += -gpu=ptxinfo,cc70,maxregcount:
endif
ifeq ($(use_device_arch),ampere)
  CXXFLAGS += -gencode=arch=compute_80,code=sm_80
  LDFLAGS  += -gencode=arch=compute_80,code=sm_80
endif
ifeq ($(use_device_arch),hopper)
  CXXFLAGS += -gencode=arch=compute_90,code=sm_90
  LDFLAGS  += -gencode=arch=compute_90,code=sm_90
endif
ifeq ($(use_device_arch),rumors)
  CXXFLAGS += -gencode=arch=compute_100,code=sm_100
  LDFLAGS  += -gencode=arch=compute_100,code=sm_100
endif

ifeq ($(use_device_arch),auto)
  # Use nvidia-smi to detect the GPU architecture
  GPU_ARCH := $(shell nvidia-smi --query-gpu=compute_cap --format=csv,noheader,nounits | awk -F. '{print "cc" $$1 $$2}')
  CXXFLAGS += -gpu=ptxinfo,$(GPU_ARCH),maxregcount:
  LDFLAGS  += -gpu=ptxinfo,$(GPU_ARCH),maxregcount:
endif

endif # use_alt_accel

CXXDEP = $(CXX)
CXXDEPFLAGS = -MM

AR = ar r
RANLIB = ranlib

TARGET_MACHINE = CUDA


#####  Fugaku_CLANG (using MPI)  ################################

else ifeq ($(strip $(target)),Fugaku_CLANG)

CPP_DEFS += -DFUJITSU_FX

ifeq ($(use_mpi),yes)
CXX = mpiFCCpx
#  CXX = mpiFCC
else ifeq ($(use_mpi),fjmpi)
CXX = mpiFCCpx
#  CXX = mpiFCC
else
CXX = FCCpx
#  CXX = FCC
endif

CXXFLAGS = -Nclang

ifeq ($(use_thread),omp)
  CXXFLAGS += -fopenmp
endif

ifeq ($(use_opt_level),low)
  OPTFLAGS = -O0
endif
ifeq ($(use_opt_level),std)
  OPTFLAGS = -O2
endif
ifeq ($(use_opt_level),high)
  OPTFLAGS = -Ofast -mllvm -inline-threshold=1000
endif

ifeq ($(use_debug),yes)
  CXXFLAGS += -g -Knolargepage
  CPP_DEFS += -DDEBUG
else
  CPP_DEFS += -DNDEBUG
endif

ifeq ($(use_qxs_arch),acle)
  CPP_DEFS += -DUSE_QXS_ACLE
endif

## NB. Fugaku_CLANG employs c++14 as default
#ifeq ($(use_cpp11),yes)
#  CXXFLAGS += -std=c++11
#endif

CXXFLAGS += $(OPTFLAGS)

LD       = $(CXX)
LDFLAGS  = $(CXXFLAGS)
LDLIBS   = -lm

ifeq ($(use_qws_library),yes)
  LDLIBS += -ltofucom -L/opt/FJSVtcs/pwrm/aarch64/lib64 -lpwr
endif

CXXDEP = $(CXX)
CXXDEPFLAGS = -MM

AR = ar r
RANLIB = ranlib

TARGET_MACHINE = Fugaku_FCC

#####  NEC SX-AURORA @KEK  #############################################

else ifeq ($(strip $(target)),SX_AURORA)

CPP_DEFS += -DNECSX
CPP_DEFS += -DUSE_STRDUP
CPP_DEFS += -DNOEXCEPTION

ifeq ($(use_mpi),yes)
  CXX = mpinc++
else
  CXX = nc++
endif

ifeq ($(use_thread),omp)
  CXXFLAGS += -fopenmp
  LDFLAGS  += -fopenmp
endif

ifeq ($(use_opt_level),low)
  OPTFLAGS = -O0 -report-all
endif
ifeq ($(use_opt_level),std)
  OPTFLAGS = -O2 -report-all
endif
ifeq ($(use_opt_level),high)
  OPTFLAGS = -O4 -report-all
endif
ifeq ($(use_opt_level),trace)
  OPTFLAGS = -O4
  CXXFLAGS += -report-all -ftrace -proginf
  LDFLAGS  += -ftrace
endif

ifeq ($(use_static),yes)
  LDFLAGS  += -static
endif

ifeq ($(use_debug),yes)
  CXXFLAGS += -g -traceback
  CPP_DEFS += -DDEBUG
else
  CPP_DEFS += -DNDEBUG
endif

ifeq ($(use_cpp11),yes)
  CXXFLAGS +=
endif

# no exceptions
CXXFLAGS += -fno-exceptions

# use restrict keyword
CXXFLAGS += -frestrict


#CXXFLAGS += -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64                      
CXXFLAGS += $(OPTFLAGS)

LD = $(CXX)
LDLIBS = -lm

CXXDEP = $(CXX)
CXXDEPFLAGS = -M

AR = ar r
RANLIB = ranlib

TARGET_MACHINE = SX_AURORA

#####
else

$(error No such a target: $(target))

endif
