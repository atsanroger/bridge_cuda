.SUFFIXES: .cpp .h .o

bridge_install_path ?= ../build
test_dir ?= ../build/tests

include $(bridge_install_path)/make.inc

## available tests
tests_avail =
tests_avail += Eigensolver
tests_avail += FFT
tests_avail += Gauge
tests_avail += GradientFlow
tests_avail += HMC
tests_avail += HotStart
tests_avail += IO
tests_avail += Mult
tests_avail += PolyakovLoop
tests_avail += QuarkNumberSusceptibility
tests_avail += RandomNumbers
tests_avail += Rational
tests_avail += SF_fAfP
tests_avail += Solver
tests_avail += Spectrum
tests_avail += TopologicalCharge
tests_avail += WilsonLoop

## main code
src_root = src

build_library_base_dir = $(bridge_install_path)
build_test_base_dir = $(test_dir)

## application codes
test=all

## check if more than one test are specified
ifneq ($(words $(test)),1)
  $(error Error: "test" should be all, none, or one of test types)
endif

test_ = $(strip $(test))

testlist=
ifeq ($(test_),all)
 testlist = $(tests_avail)
else ifeq ($(test_),none)
 testlist =
else
 ifneq ($(test_),$(filter $(test_),$(tests_avail)))
  $(error Error: unknown test $(test_))
 endif

testlist = $(test_)
endif

#--
testlist2 = $(subst QuarkNumberSusceptibility,QUARKNUMSUSCEPT,$(testlist))

test_src_dirs  =
test_src_dirs += $(src_root)
test_src_dirs += $(addprefix $(src_root)/, $(testlist))

enable_tests =
enable_tests += $(addprefix -DUSE_TEST_, $(shell echo $(testlist2) | tr a-z A-Z))

CPPFLAGS_TEST ?=
CPPFLAGS_TEST += $(test_flags)
CPPFLAGS_TEST += -I$(src_root) $(enable_tests)


program = $(build_test_base_dir)/bridge.$(test_).elf

## input files and parameters

parameter_files_dir = ParameterFiles
configuration_dir = Configurations

sample_conf_file = conf_04040408.txt

parameter_files  = $(parameter_files_dir)/main.yaml
parameter_files += $(wildcard $(addprefix $(parameter_files_dir)/,$(addsuffix /test_*.yaml, $(testlist))))

data =
data += $(wildcard $(parameter_files_dir)/parameter_rational*)



### User modification part ############################# END ###


### System part ###################################### START ###

build_test_object_dir = $(build_test_base_dir)


build_test_dirs = $(addprefix $(build_test_object_dir)/, $(test_src_dirs))

bridge_test_sources = $(wildcard $(addsuffix /*.cpp, $(test_src_dirs)))
bridge_test_objects = $(bridge_test_sources:%.cpp=$(build_test_object_dir)/%.o)
bridge_test_depends = $(addsuffix .d, $(bridge_test_objects))

depends = $(bridge_test_depends)

##--

.PHONY: all msg program cpdata all-done force

## make all

all: msg program all-done

## show summary of build options
msg:
	@( \
	echo "----------------------------------------------------------------" ; \
	echo "Target:" ; \
	echo "  test_dir      =" $(test_dir) ; \
	echo "  test          =" $(test_) ; \
	echo "  test_list     =" $(testlist) ; \
	echo "  program       =" $(program) ; \
	echo "  CPPFLAGS_TEST =" $(CPPFLAGS_TEST) ; \
	echo "----------------------------------------------------------------" ; \
	)

## build application
program: $(program) cpdata ;

$(program): $(build_test_dirs) $(bridge_test_objects)
	$(CXX) $(LDFLAGS) -o $@ $(bridge_test_objects) $(LDLIBS)
	@if [ $(test_) = 'all' ]; then cd $(build_test_base_dir) && ln -sf bridge.all.elf bridge.elf; fi

## program code
##  make directory.
$(build_test_dirs):
	@if [ ! -d $@ ]; then mkdir -p $@; echo "create directory: $@"; fi

##  find dependency and compile source file
$(build_test_object_dir)/%.o: %.cpp
	$(CXXDEP) $(CXXDEPFLAGS) $(CXXFLAGS) $(CPPFLAGS_TEST) $< \
	| sed 's|^\([^\.]*\)\.o|$(@D)/\1\.o|' > $@.d
	$(CXX) $(CXXFLAGS) $(CPPFLAGS_TEST) -c -o $@ $<

##  special cases
$(build_test_object_dir)/$(src_root)/run_test.o: $(src_root)/run_test.cpp force
	$(CXXDEP) $(CXXDEPFLAGS) $(CXXFLAGS) $(CPPFLAGS_TEST) $(enable_tests) $< \
	| sed 's|^\([^\.]*\)\.o|$(@D)/\1\.o|' > $@.d
	$(CXX) $(CXXFLAGS) $(CPPFLAGS_TEST) $(enable_tests) -c -o $@ $<

force: ;

-include $(depends)

## unzip sample conf. file and copy files
.PHONY: cpdata
cpdata:
	@echo "----------------------------------------------------------------"
	@if [ -e $(configuration_dir)/$(sample_conf_file).gz ]; then \
	  if [ ! -e $(build_test_base_dir)/$(sample_conf_file) ]; then \
	    gunzip -c $(configuration_dir)/$(sample_conf_file).gz > \
	      $(build_test_base_dir)/$(sample_conf_file) ; \
	    echo copy $(configuration_dir)/$(sample_conf_file); \
	  else \
	    echo skip $(configuration_dir)/$(sample_conf_file); \
	  fi; \
	fi;
	@for f in $(data); do \
	  f2=`basename $$f`; \
	  if [ ! -e $(build_test_base_dir)/$$f2 ]; then \
	    echo copy $$f; \
	    cp $$f $(build_test_base_dir); \
	  elif [ $$f -nt $(build_test_base_dir)/$$f2 ]; then \
	    echo update $$f; \
	    mv $(build_test_base_dir)/$$f2 $(build_test_base_dir)/$$f2.bak; \
	    cp $$f $(build_test_base_dir); \
	  else \
	    echo skip $$f; \
	  fi; \
	done
	@for f in $(parameter_files); do \
	  f2=`basename $$f`; \
	  if [ ! -e $(build_test_base_dir)/$$f2 ]; then \
	    echo copy $$f; \
	    cp $$f $(build_test_base_dir); \
	  elif [ $$f -nt $(build_test_base_dir)/$$f2 ]; then \
	    echo update $$f; \
	    mv $(build_test_base_dir)/$$f2 $(build_test_base_dir)/$$f2.bak; \
	    cp $$f $(build_test_base_dir); \
	  else \
	    echo skip $$f; \
	  fi; \
	done
	@echo "----------------------------------------------------------------"

all-done: $(program)
	@echo ""
	@echo "$(program) is successfully compiled."
	@echo ""

## exec. test run (single ver. only)
.PHONY: run
run:
	@if [ -e $(program) ]; then \
	  echo "run program ..."; \
	  ( \
	    cd $(build_test_base_dir); \
	    program_=`basename $(program)`; \
	    ./$$program_ -a \
	  ) \
	fi

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
	if [ -d "$(build_test_object_dir)/$(src_root)" ]; then rm -rf $(build_test_object_dir)/$(src_root); fi
	rm -f $(program)
	if [ $(test_) = 'all' ]; then cd $(build_test_base_dir); rm -f bridge.elf; fi

#	
#	

################################################################################
#        @author  Shinji Motoki (smotoki) 
#                 $LastChangedBy$
#        @date    $LastChangedDate::                      $
#        @version $LastChangedRevision$
################################################################################
