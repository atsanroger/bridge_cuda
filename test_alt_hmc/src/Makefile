.SUFFIXES: .cpp .h .o

BRIDGE_INSTALL_PATH = ../../build

include $(BRIDGE_INSTALL_PATH)/make.inc

CXXFLAGS += -I./

sources = main.cpp \
 Tests/test.cpp \
 Tests/test_HMC_Staggered_Leapfrog_alt.cpp \
 Tests/test_HMC_Domainwall_Nf2_PV.cpp \
 Tests/test_HMC_Domainwall_Nf2_PV_nosmr.cpp \
 Tests/test_HMC_Domainwall_Nf2_PV_alt_Field.cpp \
 Tests/test_HMC_Domainwall_Nf2_PV_alt_nosmr.cpp \
 Tests/test_HMC_Domainwall_Nf2_PV_alt.cpp \
 Tests/test_HMC_Domainwall_Nf2_PV_alt_eo_nosmr.cpp \
 Tests/test_HMC_Domainwall_Nf2_PV_alt_eo.cpp \
 Tests/test_HMC_Domainwall_Nf2p1_PV_alt.cpp \
 Tests/test_HMC_Domainwall_Nf2p1_PV_alt_eo.cpp \
 Tests/test_HMC_Wilson_Nf2.cpp \
 Tests/test_HMC_Wilson_Nf2_nosmr.cpp \
 Tests/test_HMC_Wilson_Nf2_alt_Field.cpp \
 Tests/test_HMC_Wilson_Nf2_alt_nosmr.cpp \
 Tests/test_HMC_Wilson_Nf2_alt.cpp \
 Tests/test_HMC_Wilson_Nf2_alt_eo.cpp \
 Tests/test_HMC_Wilson_Nf2_alt_eo_ratio.cpp \
 Tests/test_HMC_Wilson_Nf2_alt_eo_nosmr.cpp \
 Tests/test_HMC_Wilson_RHMC_Nf2p1.cpp \
 Tests/test_HMC_Wilson_RHMC_Nf2p1_alt.cpp \
 Tests/test_HMC_Wilson_RHMC_Nf2p1_alt_eo.cpp \
 Tests/test_HMC_Wilson_RHMC_Nf2p1_alt_eo_ratio.cpp \
 Tests/job_Utils.cpp


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

