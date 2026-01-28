from pybridge import *

setup(lattice_size=[4,4,4,8], grid_size=[], verbose_level="Detailed")

nvol=Common.Nvol
ndim=Common.Ndim

u=Field_G(nvol, ndim)

#--- Load gauge config

gconf=GaugeConfig("Text")
gconf.read_file(u, "conf_04040408.txt")

# // shift

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

smear=Director_Smear(Smear("APE",
                           Projection("Stout_SU3",
                                      maximum_number_of_iteration=1000),
                           rho_uniform=0.1),
                     number_of_smearing=0)

#--- Propagator calculator
#
def calc_prop(kappa, slope):
    fopr=Fopr_Smeared(Fopr("Clover",
                           gamma_matrix_type="Dirac",
                           hopping_parameter=kappa,
                           clover_coefficient=1.0,
                           boundary_condition=[1,1,1,-1]),
                      smear)

    solver=Solver("BiCGStab_Cmplx",
                  fopr,
                  maximum_number_of_iteration=1000,
                  maximum_number_of_restart=40,
                  convergence_criterion_squared=1.0e-28,
                  use_initial_guess="false",
                  Omega_tolerance=0.7,
                  verbose_level="General")

    fprop=Fprop_Standard_lex(solver)

    source=Source("Exponential",
                  source_position=[0,0,0,0],
                  slope=slope,
                  power=2.0)

    fprop.set_config(u)

    nc=Common.Nc
    nd=Common.Nd

    pp=[0] * (nc*nd)

    for ispin in range(nd):
        for icolor in range(nc):
            icd=icolor + nc * ispin
        
            b=Field_F(nvol)
            source.set(b, icolor=icolor, ispin=ispin)
    
            x=b.clone()

            nconv, diff = fprop.invert_D(x, b)
            print("ic= {}, id= {}, nconv= {}, diff= {}".format(icolor, ispin, nconv, diff))

            pp[icd] = x

    print(pp)

    return pp


#--- calculate propagators

pp1=calc_prop(kappa=0.12, slope=0.25)
pp2=calc_prop(kappa=0.115, slope=0.5)

#--- 2pt correlator

corr=Corr2pt_4spinor(GammaMatrixSet("Dirac"))

print("=== quark1 - quark1")
corr.meson_all(pp1, pp1)

print("=== quark2 - quark2")
corr.meson_all(pp2, pp2)

print("=== quark1 - quark2")
corr.meson_all(pp1, pp2)
