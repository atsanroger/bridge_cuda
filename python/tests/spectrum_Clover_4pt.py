from pybridge import *

setup(lattice_size=[4,4,4,8], grid_size=[],
      random_number_type="Mseries",
      random_number_seed=1234567,
      verbose_level="Detailed")

nvol=Common.Nvol
ndim=Common.Ndim

u=Field_G(nvol, ndim)

#--- Load gauge config

gconf=GaugeConfig("Text")
gconf.read_file(u, "conf_04040408.txt")

#--- Shift origin

shift_origin=[0,0,0,0]

shift=ShiftField_lex()

for idir in range(len(shift_origin)):
    if shift_origin[idir] == 0:
        pass
    elif shift_origin[idir] > 0:
        for i in range(shift_origin[idir]):
            ut=u.clone()
            shift.backward(ut,u,idir)
            copy(u,ut)
    else:
        for i in range(-shift_origin[idir]):
            ut=u.clone()
            shift.forward(ut,u,idir)
            copy(u,ut)

#--- Gauge Fixing

gfix=GaugeFixing("Coulomb",
                 maximum_number_of_iteration=5000,
                 number_of_naive_iteration=50,
                 interval_of_measurement=10,
                 iteration_to_reset=1000,
                 convergence_criterion_squared=1.0e-24,
                 overrelaxation_parameter=1.6,
                 verbose_level="Detailed")

u2=u.clone()
gfix.fix(u2,u)
u=u2

#--- Smearing

smear=Director_Smear(Smear("APE",
                           Projection("Stout_SU3"),
                           rho_uniform=0.1),
                     number_of_smearing=0)

#--- Propagator calculator
#
def calc_prop(fopr, solver, source):

    fprop=Fprop_Standard_lex(solver)
    fprop.set_config(u)

    nc=Common.Nc
    nd=Common.Nd

    print("Solving quark propagator:")
    print("  color spin   Nconv      diff           diff2")

    pp=[0] * (nc*nd)

    for ispin in range(nd):
        for icolor in range(nc):
            icd=icolor + nc * ispin
        
            b=Field_F(nvol)
            source.set(b, icolor=icolor, ispin=ispin)
    
            x=b.clone()
            nconv, diff = fprop.invert_D(x, b)

            y=b.clone()
            copy(y,b)

            fopr.set_mode("D")
            fopr.mult(y,x)
            axpy(y, -1.0, b)
            diff2=y.norm2()/b.norm2()

            print("   {:2d}   {:2d}   {:6d}   {:12.4e}   {:12.4e}"
                  .format(icolor, ispin, nconv, diff, diff2))

            pp[icd] = x

    print(pp)

    return pp


#--- calculate propagators

fopr1 = Fopr_Smeared(Fopr("Clover",
                          gamma_matrix_type="Dirac",
                          hopping_parameter=0.12,
                          clover_coefficient=1.0,
                          boundary_condition=[1,1,1,-1]),
                     smear)

solver1 = Solver("BiCGStab_Cmplx",
                 fopr1,
                 maximum_number_of_iteration=100,
                 maximum_number_of_restart=40,
                 convergence_criterion_squared=1.0e-28,
                 use_initial_guess="false",
                 Omega_tolerance=0.70,
                 verbose_level="General")

source1 = Source("Random",
                 source_position=[0,0,0,0],
                 source_momentum=[1,1,1,0],
                 noise_type="Z2")

pp1 = calc_prop(fopr1, solver1, source1)



fopr2 = Fopr_Smeared(Fopr("Clover",
                          gamma_matrix_type="Dirac",
                          hopping_parameter=0.12,
                          clover_coefficient=1.0,
                          boundary_condition=[1,1,1,-1]),
                     smear)

solver2 = Solver("BiCGStab_Cmplx",
                 fopr2,
                 maximum_number_of_iteration=100,
                 maximum_number_of_restart=40,
                 convergence_criterion_squared=1.0e-28,
                 use_initial_guess="false",
                 Omega_tolerance=0.70,
                 verbose_level="General")

source2 = Source("Random",
                 source_position=[0,0,0,0],
                 source_momentum=[1,1,1,0],
                 noise_type="Z2")

pp2 = calc_prop(fopr2, solver2, source2)


#--- 4pt correlator

corr_4pt = Corr4pt_4spinor(GammaMatrixSet("Dirac"))

print("==== 4-point correlator ====")
corr_4pt.meson_all(pp1, pp1, pp2, pp2)


#--- 2pt correlator

corr=Corr2pt_4spinor(GammaMatrixSet("Dirac"))

print("=== quark1 - quark1")
corr.meson_all(pp1, pp1)

print("=== quark2 - quark2")
corr.meson_all(pp2, pp2)

print("=== quark1 - quark2")
corr.meson_all(pp1, pp2)
