from pybridge import *

setup(lattice_size=[4,4,4,8], grid_size=[], verbose_level="General")

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

#--- Propagator calculation
#
fopr_w=Fopr_Wilson(gamma_matrix_type="Dirac",
                              hopping_parameter=0.5/(4.0-1.6),
                              boundary_condition=[1,1,1,1])
fopr=Fopr_Overlap(fopr_w,
                  quark_mass=0.2,
                  domain_wall_height=1.6,
                  number_of_poles=16,
                  lower_bound=0.01,
                  upper_bound=8.0,
                  maximum_number_of_iteration=1000,
                  convergence_criterion_squared=1.0e-20,
                  boundary_condition=[1,1,1,1])

solver=Solver("CGNR",
              fopr,
              maximum_number_of_iteration=100,
              maximum_number_of_restart=40,
              convergence_criterion_squared=1.0e-20,
              use_initial_guess="false"
)

fprop=Fprop_Standard_lex(solver)

fprop.set_config(u)

source=Source("Local",
              source_position=[0,0,0,0]
)


nc=Common.Nc
nd=Common.Nd

pp=[0] * (nc*nd)

print("Solving quark propagator:")
print("  color spin   Nconv      diff           diff2")

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
        diff2=y.norm2()
        
        print("   {:2d}   {:2d}   {:6d}   {:12.4e}   {:12.4e}"
              .format(icolor, ispin, nconv, diff, diff2))

        pp[icd] = x

print(pp)

#--- 2pt correlator

corr=Corr2pt_4spinor(GammaMatrixSet("Dirac"))

corr.meson_all(pp, pp)
