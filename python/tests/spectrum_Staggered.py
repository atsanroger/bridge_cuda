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

#--- Split in even-odd

ueo=Field_G(nvol, ndim)

Index_eo().convert_field(ueo, u)

print("m_conf norm2 = {:20.16e}".format(u.norm2()))
print("Ueo    norm2 = {:20.16e}".format(ueo.norm2()))

#--- Build objects
#
fopr=Fopr_Staggered_eo(quark_mass=0.1,
                       boundary_condition=[1,1,1,-1])

fopr.set_config(ueo)

fprop=Fprop_Staggered_eo(fopr,
                         solver_type="GMRES_m_Cmplx",
                         maximum_number_of_iteration=100,
                         maximum_number_of_restart=40,
                         convergence_criterion_squared=1.0e-28,
                         use_initial_guess="false",
                         number_of_orthonormal_vectors=2)

source=Source_Staggered_Wall(source_position=0)



#--- Propagator calculation
#
nc=Common.Nc

pp=[0] * (nc*2)

for i_src in range(2):
    for ic in range(nc):
        print("id = {}, i_src = {}".format(ic, i_src))

        icd=ic + nc * i_src
        
        b=Field_F_1spinor(nvol)
        source.set(b, icolor=ic, isource=i_src)
    
        x=Field_F_1spinor(nvol)
        fprop.invert_D(x, b)

        pp[icd] = x

print(pp)


#--- 2pt correlator
#
corr=Corr2pt_Staggered()

r=corr.meson(pp, pp)

print("PS <-- PS correlator:")

for t in range(len(r)):
    print("  {:4d} {:20.12e}".format(t, r[t]))


