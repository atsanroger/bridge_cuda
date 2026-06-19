.SUFFIXES: .py .log

PYTHON=env PYTHONPATH=../src python3
#PYTHON=env PYTHONPATH=../src mpiexec -np 4 python3

.PHONY: Eigensolver FFT Gauge GradientFlow HMC HotStart IO Mult PolyakovLoop QuarkNumberSusceptibility RandomNumbers Rational SF_fAfP Solver Spectrum TopologicalCharge WilsonLoop

TESTS= \
	Eigensolver \
	FFT \
	Gauge \
	GradientFlow \
	HMC \
	HotStart \
	IO \
	Mult \
	PolyakovLoop \
	QuarkNumberSusceptibility \
	RandomNumbers \
	Rational \
	SF_fAfP \
	Solver \
	Spectrum \
	TopologicalCharge \
	WilsonLoop


all: $(TESTS)

Eigensolver: \
	eigensolver.log \
	eigensolver_chebyshev.log

FFT: \
	fft_check.log

Gauge: \
	gauge_plaquette.log \
	gauge_shift.log \
	gauge_energy_density.log

GradientFlow: \
	gradientflow_nf0_rk1.log \
	gradientflow_nf0_rk2.log \
	gradientflow_nf0_rk3.log \
	gradientflow_nf0_rk4.log \
	gradientflow_nf0_rk3_adapt.log \
	gradientflow_emt_gauge.log \
	gradientflow_emt_fermion1pt.log \
	gradientflow_emt_fermion2pt.log

gradientflow_nf0_rk1.log: gradientflow_nf0.py
	$(PYTHON) $< -v --rk=1 | tee $@

gradientflow_nf0_rk2.log: gradientflow_nf0.py
	$(PYTHON) $< -v --rk=2 | tee $@

gradientflow_nf0_rk3.log: gradientflow_nf0.py
	$(PYTHON) $< -v --rk=3 | tee $@

gradientflow_nf0_rk4.log: gradientflow_nf0.py
	$(PYTHON) $< -v --rk=4 | tee $@

gradientflow_nf0_rk3_adapt.log: gradientflow_nf0_adapt.py
	$(PYTHON) $< | tee $@

HMC: \
	hmc_clover_nf2.log \
	hmc_clover_nf2_hyp.log \
	hmc_clover_nf2p1.log \
	hmc_wilson_nf2.log \
	hmc_wilson_nf2_leapfrog.log \
	hmc_quenched.log \
	hmc_quenched_leapfrog.log \
	hmc_clover_isochemical_nf2.log \
	hmc_clover_isochemical_nf2p1.log \
	hmc_clover_nf2_eo.log \
	hmc_clover_nf2p1_eo.log \
	hmc_staggered_nf2_eo.log \
	hmc_domainwall_nf2.log \
	hmc_domainwall_nf2_PV.log \
	hmc_overlap_nf2.log \
	hmc_wilson_twistedmass_nf2.log \
	hmc_wilson_nf2_topologyfix.log

HotStart: \
	hotstart.log \
	hotstart_det.log

IO: \
	io_data_text.log \
	io_gaugeconfig.log

Mult: \
	mult.log \
	mult_eo.log \
	mult_domainwall.log \
	mult_overlap.log \
	mult_wilson_twistedmass.log

PolyakovLoop: \
	polyakov_loop.log

QuarkNumberSusceptibility: \
	quarknum.log

RandomNumbers: \
	rand_field.log \
	rand_uniform.log

Rational: \
	rational_approx.log \
	rational_inverse.log \
	rational_smeared.log 

SF_fAfP:

Solver: \
	solver.log \
	solver_shiftcg.log

Spectrum: \
	spectrum_Wilson.log \
	spectrum_Wilson_shift_origin.log \
	spectrum_WilsonGeneral.log \
	spectrum_Clover.log \
	spectrum_CloverGeneral.log \
	spectrum_Clover_4pt.log \
	spectrum_Clover_fieldio.log \
	spectrum_Clover_eo.log \
	spectrum_CRS.log \
	spectrum_Staggered.log \
	spectrum_Domainwall.log \
	spectrum_Overlap.log \
	spectrum_Wilson_TwistedMass.log

TopologicalCharge: \
	topological_charge.log

WilsonLoop: \
	wilsonloop.log


.py.log:
	$(PYTHON) $< | tee $@

clean:
	@echo "unimplemented."
