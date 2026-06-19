.SUFFIXES: .cpp .h .o

BRIDGE_INSTALL_PATH = ../../build

include $(BRIDGE_INSTALL_PATH)/make.inc

CXXFLAGS += -I./

sources = main.cpp \
 Tests/job_Utils.cpp \
 Tests/test.cpp \
 Tests/spectrum_alt.cpp \
 Tests/test_alt_Corelib.cpp \
 Tests/test_alt_SIMD2.cpp \
 Tests/test_alt_Vector.cpp \
 Tests/test_alt_QXS.cpp \
 Tests/test_alt_Accel.cpp


objects = $(sources:%.cpp=%.o)

program = bridge.elf

all: msg $(program)

$(program): $(objects)
	$(LD) $(LDFLAGS) -o $@ $(objects) $(LDLIBS)

%.o: %.cpp $(headers)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

msg:
	@cat $(BRIDGE_INSTALL_PATH)/message.txt

clean:
	rm -f $(program) $(objects)

