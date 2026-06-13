.SUFFIXES: .cpp .h .o

#POSTFIX = _mpi
#POSTFIX = _mpi_debug
#POSTFIX = _low_debug
#POSTFIX = _low
#POSTFIX = _high
#POSTFIX = _mpi_low_debug
#POSTFIX = _mpi_high
#POSTFIX = _debug
POSTFIX =

BRIDGE_INSTALL_PATH = ../build$(POSTFIX)

#lib_alt_flag = -lbridge_alt

include $(BRIDGE_INSTALL_PATH)/make.inc

# for alternative: no need anymore
#alt_include_dirs = \
# $(BRIDGE_INSTALL_PATH)/include/bridge/lib_alt\
# $(BRIDGE_INSTALL_PATH)/include/bridge/lib_alt_Simd\

# -I$(BRIDGE_INSTALL_PATH)/include/bridge/lib_alt/Fopr\
# -I$(BRIDGE_INSTALL_PATH)/include/bridge/lib_alt_Simd/Field\
# -I$(BRIDGE_INSTALL_PATH)/include/bridge/lib_alt_Simd/Fopr\


#DEBUGFLAGS = -g -DDEBUG

#headers = test.h

src_dirs = ./
build_base_dir = build$(POSTFIX)
build_obj_dir = $(build_base_dir)/obj

CXXFLAGS += $(addprefix -I, $(src_dirs)) $(addprefix -I, $(alt_include_dirs))
CXXDEPFLAGS += $(BRIDGE_INCLUDE)
CXXDEPFLAGS += $(addprefix -I, $(src_dirs)) $(addprefix -I, $(alt_include_dirs)) -Isrc/Tests

#DEBUGFLAGS = -g -DDEBUG


sources = $(wildcard $(addsuffix /*.cpp, $(src_dirs)))
objects = $(sources:%.cpp=$(build_obj_dir)/%.o)
depends = $(addsuffix .d, $(objects))
build_dirs = $(addprefix $(build_obj_dir)/, $(src_dirs))

# object file to show version (always compiled)
src_show_version =show_version.cpp
obj_show_version =$(build_obj_dir)/show_version.o

## version info
#GIT_COMMIT_AUTHOR := $(shell git log -1 --format='Author: %an <%ae>')
#GIT_COMMIT_ID := $(shell git log -1 --format='commit %H')
#GIT_COMMIT_DATE := $(shell git log -1 --format='%ad')
BRIDGE_VERSION := $(shell cat ../VERSION)

#VERSION_INFO= -D_BRIDGE_VERSION="$(BRIDGE_VERSION)" -D_COMMIT_AUTHOR="$(GIT_COMMIT_AUTHOR)"  -D_COMMIT_ID="$(GIT_COMMIT_ID)"  -D_COMMIT_DATE="$(GIT_COMMIT_DATE)"

VERSION_INFO= -D_BRIDGE_VERSION="$(BRIDGE_VERSION)"


program = $(build_base_dir)/bridge_test_alt

all: $(program) msg

$(program): $(build_dirs) $(objects) $(obj_show_version)
	$(LD) $(LDFLAGS) -o $@ $(objects) $(LDLIBS)

## find dependency and compile source file
$(build_obj_dir)/%.o: %.cpp
	$(CXXDEP) $(CPPFLAGS) $(CXXDEPFLAGS) $< \
	| sed 's|^\([^\.]*\)\.o|$(@D)/\1\.o|' > $@.d
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(DEBUGFLAGS) -c -o $@ $<

.PHONY: $(src_show_version)
$(obj_show_version): $(src_show_version)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(DEBUGFLAGS) $(VERSION_INFO) -c -o $@ $<

msg:
	@cat $(BRIDGE_INSTALL_PATH)/message.txt
	@echo $(build_dirs)

lib_alt:
	@cd ../ ; make lib_alt

## make directory.
$(build_dirs):
	@if [ ! -d $@ ]; then mkdir -p $@; echo "create directory: $@"; fi


clean:
	rm -f $(program) $(objects)

-include $(depends)
