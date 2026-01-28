.SUFFIXES: .cpp .h .hh .o

build_dir ?= ./build

## main code
src_root = src

build_library_base_dir = $(build_dir)

## library code

lib_bridge = $(build_library_base_dir)/libbridge.a
lib_bridge_include = $(build_library_base_dir)/include/bridge

#-- directories containing source files of bridge library
lib_bridge_src_root = $(src_root)/lib

lib_bridge_src_dirs =
lib_bridge_src_dirs += $(lib_bridge_src_root)
lib_bridge_src_dirs += $(lib_bridge_src_root)/Action
lib_bridge_src_dirs += $(lib_bridge_src_root)/Action/Fermion
lib_bridge_src_dirs += $(lib_bridge_src_root)/Action/Gauge
lib_bridge_src_dirs += $(lib_bridge_src_root)/Communicator
lib_bridge_src_dirs += $(lib_bridge_src_root)/Communicator/$(SRCCOMMUN)
lib_bridge_src_dirs += $(lib_bridge_src_root)/Eigen
lib_bridge_src_dirs += $(lib_bridge_src_root)/Field
lib_bridge_src_dirs += $(lib_bridge_src_root)/Field/$(SRCOPT)
lib_bridge_src_dirs += $(lib_bridge_src_root)/Fopr
lib_bridge_src_dirs += $(lib_bridge_src_root)/Fopr/Org
lib_bridge_src_dirs += $(lib_bridge_src_root)/Fopr/Imp
#lib_bridge_src_dirs += $(lib_bridge_src_root)/Fopr/$(SRCOPT2)
lib_bridge_src_dirs += $(lib_bridge_src_root)/Force/Fermion
lib_bridge_src_dirs += $(lib_bridge_src_root)/Force/Fermion/$(SRCOPT)
lib_bridge_src_dirs += $(lib_bridge_src_root)/Force/Gauge
lib_bridge_src_dirs += $(lib_bridge_src_root)/HMC
lib_bridge_src_dirs += $(lib_bridge_src_root)/IO
lib_bridge_src_dirs += $(lib_bridge_src_root)/Measurements/Fermion
lib_bridge_src_dirs += $(lib_bridge_src_root)/Measurements/Fermion/$(SRCOPT)
lib_bridge_src_dirs += $(lib_bridge_src_root)/Measurements/Gauge
lib_bridge_src_dirs += $(lib_bridge_src_root)/Parameters
lib_bridge_src_dirs += $(lib_bridge_src_root)/ResourceManager
lib_bridge_src_dirs += $(lib_bridge_src_root)/ResourceManager/$(SRCTHREAD)
lib_bridge_src_dirs += $(lib_bridge_src_root)/Smear
lib_bridge_src_dirs += $(lib_bridge_src_root)/Solver
lib_bridge_src_dirs += $(lib_bridge_src_root)/Tools

ifeq ($(use_alternative),yes)
 lib_bridge_alt_root = $(src_root)/lib_alt
 alt_dirs += $(lib_bridge_alt_root)
 alt_dirs += $(lib_bridge_alt_root)/Action/Fermion
 alt_dirs += $(lib_bridge_alt_root)/Eigen
 alt_dirs += $(lib_bridge_alt_root)/Field
 alt_dirs += $(lib_bridge_alt_root)/Fopr
 alt_dirs += $(lib_bridge_alt_root)/Force/Fermion
 alt_dirs += $(lib_bridge_alt_root)/Measurements/Gauge
 alt_dirs += $(lib_bridge_alt_root)/Measurements/Fermion
 alt_dirs += $(lib_bridge_alt_root)/Solver
 alt_dirs += $(lib_bridge_alt_root)/Smear
 CPPFLAGS += -I$(lib_bridge_alt_root)


ifeq ($(use_alt_qxs),yes)
 lib_bridge_alt_qxs_root = $(src_root)/lib_alt_QXS
 alt_dirs += $(lib_bridge_alt_qxs_root)
 alt_dirs += $(lib_bridge_alt_qxs_root)/inline
 alt_dirs += $(lib_bridge_alt_qxs_root)/extra
 alt_dirs += $(lib_bridge_alt_qxs_root)/ResourceManager
 alt_dirs += $(lib_bridge_alt_qxs_root)/Action/Fermion
 alt_dirs += $(lib_bridge_alt_qxs_root)/Eigen
 alt_dirs += $(lib_bridge_alt_qxs_root)/Field
 alt_dirs += $(lib_bridge_alt_qxs_root)/Fopr
 alt_dirs += $(lib_bridge_alt_qxs_root)/Force/Fermion
 alt_dirs += $(lib_bridge_alt_qxs_root)/Measurements/Gauge
 alt_dirs += $(lib_bridge_alt_qxs_root)/Measurements/Fermion
 alt_dirs += $(lib_bridge_alt_qxs_root)/Solver
 alt_dirs += $(lib_bridge_alt_qxs_root)/Smear
 alt_dirs += $(lib_bridge_alt_qxs_root)/BridgeQXS
 CPPFLAGS += -I$(lib_bridge_alt_qxs_root)
endif

ifeq ($(use_alt_accel),yes)
 lib_bridge_alt_accel_root = $(src_root)/lib_alt_Accel
 alt_dirs += $(lib_bridge_alt_accel_root)
 alt_dirs += $(lib_bridge_alt_accel_root)/inline
 alt_dirs += $(lib_bridge_alt_accel_root)/extra
 alt_dirs += $(lib_bridge_alt_accel_root)/ResourceManager
 alt_dirs += $(lib_bridge_alt_accel_root)/Action/Fermion
 alt_dirs += $(lib_bridge_alt_accel_root)/Eigen
 alt_dirs += $(lib_bridge_alt_accel_root)/Field
 alt_dirs += $(lib_bridge_alt_accel_root)/Fopr
 alt_dirs += $(lib_bridge_alt_accel_root)/Force/Fermion
 alt_dirs += $(lib_bridge_alt_accel_root)/Measurements/Gauge
 alt_dirs += $(lib_bridge_alt_accel_root)/Measurements/Fermion
 alt_dirs += $(lib_bridge_alt_accel_root)/Solver
 alt_dirs += $(lib_bridge_alt_accel_root)/BridgeACC
 CPPFLAGS += -I$(lib_bridge_alt_accel_root)
endif

ifeq ($(use_alt_vector),yes)
 lib_bridge_alt_vector_root = $(src_root)/lib_alt_Vector
 alt_dirs += $(lib_bridge_alt_vector_root)
 alt_dirs += $(lib_bridge_alt_vector_root)/inline
 alt_dirs += $(lib_bridge_alt_vector_root)/ResourceManager
 alt_dirs += $(lib_bridge_alt_vector_root)/Action/Fermion
 alt_dirs += $(lib_bridge_alt_vector_root)/Eigen
 alt_dirs += $(lib_bridge_alt_vector_root)/Field
 alt_dirs += $(lib_bridge_alt_vector_root)/Fopr
 alt_dirs += $(lib_bridge_alt_vector_root)/Force/Fermion
 alt_dirs += $(lib_bridge_alt_vector_root)/Measurements/Gauge
 alt_dirs += $(lib_bridge_alt_vector_root)/Measurements/Fermion
 alt_dirs += $(lib_bridge_alt_vector_root)/Solver
 alt_dirs += $(lib_bridge_alt_vector_root)/Smear
 alt_dirs += $(lib_bridge_alt_vector_root)/BridgeVec
 CPPFLAGS += -I$(lib_bridge_alt_vector_root)
endif

 lib_bridge_src_dirs += $(alt_dirs)
endif

##--

#CPPFLAGS = -I$(src_root) -I$(lib_bridge_src_root)
#CPPFLAGS += $(EXTRA_INCLUDES)
#CPPFLAGS += $(CPP_DEFS)
CPPFLAGS += -Isrc -Isrc/lib

##--
build_library_dirs = $(addprefix $(build_library_base_dir)/, $(lib_bridge_src_dirs))

lib_bridge_sources = $(wildcard $(addsuffix /*.cpp, $(lib_bridge_src_dirs)))
lib_bridge_objects = $(lib_bridge_sources:%.cpp=$(build_library_base_dir)/%.o)
lib_bridge_depends = $(addsuffix .d, $(lib_bridge_objects))

#lib_bridge_flag = -L$(build_library_base_dir) -Wl,--whole-archive -lbridge -Wl,--no-whole-archive
lib_bridge_flag = -L$(build_library_base_dir) -lbridge

make_inc_file = $(build_library_base_dir)/make.inc
msg_file = $(build_library_base_dir)/message.txt

depends = $(bridge_test_depends) $(lib_bridge_depends)

##--

.PHONY: all msg all-done
.PHONY: lib cpinclude make-inc

## make all

all:
	@echo bridge++ library: run \"make\" in top directory.

## build library
lib: msg $(lib_bridge) cpinclude all-done ;

$(lib_bridge): $(build_library_dirs) $(lib_bridge_objects)
	$(AR) $@ $(lib_bridge_objects)
	$(RANLIB) $@

## show summary of build options
msg:
##  make directory.
	@( \
	echo "----------------------------------------------------------------"; \
	echo "Compilation Environment Summary" ; \
	echo ; \
	echo "Options:" ; \
	echo "  target              =" $(target) "(" $(TARGET_MACHINE) ")" ; \
	echo "  use_mpi             =" $(use_mpi) ; \
	echo "  use_thread          =" $(use_thread) ; \
	echo "  use_opt_code        =" $(use_opt_code) ; \
	echo "  use_opt_level       =" $(use_opt_level) ; \
	echo "  use_static          =" $(use_static) ; \
	echo "  use_debug           =" $(use_debug) ; \
	echo "  use_gauge_group     =" $(use_gauge_group) ; \
	echo "  use_complex         =" $(use_complex) ; \
	echo "  use_cpp11           =" $(use_cpp11) ; \
	echo "  use_testmanager     =" $(use_testmanager) ; \
	echo "  use_lime_library    =" $(use_lime_library) ; \
	echo "  use_sfmt_library    =" $(use_sfmt_library) ; \
	echo "  use_ntl_library     =" $(use_ntl_library) ; \
	echo "  use_fftw_library    =" $(use_fftw_library) ; \
	echo "  use_yamlcpp_library =" $(use_yamlcpp_library) ; \
	echo "  use_tinyxml2_library=" $(use_tinyxml2_library) ; \
	echo ; \
	echo "Compiler and options:" ; \
	echo "  CXX         =" $(CXX) ; \
	echo "  CXXFLAGS    =" $(CXXFLAGS) ; \
	echo "  LD          =" $(LD) ; \
	echo "  LDFLAGS     =" $(LDFLAGS) ; \
	echo "  LDLIBS      =" $(LDLIBS) ; \
	echo "  CXXDEP      =" $(CXXDEP) ; \
	echo "  CXXDEPFLAGS =" $(CXXDEPFLAGS) ; \
	echo "  CPPFLAGS    =" $(CPPFLAGS) ; \
	echo "----------------------------------------------------------------" ; \
	)

## library code
##  make directory.
$(build_library_dirs):
	@if [ ! -d $@ ]; then mkdir -p $@; echo "create directory: $@"; fi

##  find dependency and compile source file
$(build_library_base_dir)/%.o: %.cpp
	$(CXXDEP) $(CXXDEPFLAGS) $(CPPFLAGS) $(CXXFLAGS) $< \
	| sed 's|^\([^\.]*\)\.o|$(@D)/\1\.o|' > $@.d
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<


-include $(depends)

## unzip sample conf. file and copy files

cpinclude:
	@mkdir -p $(lib_bridge_include)
	@echo copy library header files ...
	@echo $(src_root) $(lib_bridge_src_dirs) | \
        while read i; \
        do \
	  find -L $$i -maxdepth 1 -type f \( -name '*.h' -o -name '*.hh' \) ; \
	done | \
	sed -e 's,^$(src_root)/,,' | \
	( \
	  cd $(src_root); \
	  cpio -o --quiet \
	) | ( \
	  cd $(lib_bridge_include); \
	  cpio -idmuv --quiet \
	)

all-done: $(program)
	@echo ""
	@echo "Bridge library is successfully compiled."
	@echo ""

## generate doxygen html files to ./docs/html
.PHONY: doxygen
doxygen:
	@if [ -e ./docs/doxygen.conf ]; then \
	  echo "generating doxygen documents ..."; \
	  ( \
	    cd ./docs; \
	    doxygen doxygen.conf \
	  ) \
	fi

## file clean up
.PHONY: clean
clean:
	if [ -d "$(build_library_base_dir)/$(lib_bridge_src_root)" ]; then rm -rf $(build_library_base_dir)/$(lib_bridge_src_root); fi

real-clean: clean
	rm -f $(lib_bridge) $(make_inc_file) $(msg_file)
	if [ -d "$(lib_bridge_include)" ]; then rm -rf $(lib_bridge_include); fi
#	
#	

################################################################################
#        @author  Shinji Motoki (smotoki) 
#                 $LastChangedBy$
#        @date    $LastChangedDate::                      $
#        @version $LastChangedRevision$
################################################################################
