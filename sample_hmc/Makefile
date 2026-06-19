.SUFFIXES: .cpp .h .o

BRIDGE_INSTALL_PATH = ../build

include $(BRIDGE_INSTALL_PATH)/make.inc

headers =
sources = sample_HMC.cpp
objects = $(sources:%.cpp=%.o)

program = sample_HMC.elf

all: msg $(program)

$(program): $(objects)
	$(LD) $(LDFLAGS) -o $@ $(objects) $(LDLIBS)

%.o: %.cpp $(headers)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

msg:
	@cat $(BRIDGE_INSTALL_PATH)/message.txt

clean:
	rm -f $(program) $(objects)

