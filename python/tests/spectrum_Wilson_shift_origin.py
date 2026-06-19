from pybridge import *

setup(lattice_size=[4,4,4,8], grid_size=[], verbose_level="Detailed")

nvol=Common.Nvol
ndim=Common.Ndim

u=Field_G(nvol, ndim)

#--- Load gauge config

gconf=GaugeConfig("Text")
gconf.read_file(u, "conf_04040408.txt")

#--- Shift origin

#shift_origin=[0,0,0,2]
#source_position=[0,0,0,6]
shift_origin=[0,0,0,0]
source_position=[0,0,0,0]

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

smear=Director_Smear(Smear("APE",
                           Projection("Stout_SU3"),
                           rho_uniform=0.1),
                     number_of_smearing=0)

fopr=Fopr_Smeared(Fopr("Wilson",
                       gamma_matrix_type="Chiral",
                       hopping_parameter=0.12,
                       boundary_condition=[1,1,1,-1]),
                  smear)

solver=Solver("CGNR",
              fopr,
              maximum_number_of_iteration=100,
              maximum_number_of_restart=40,
              convergence_criterion_squared=1.0e-28,
              use_initial_guess="false")

fprop=Fprop_Standard_lex(solver)

fprop.set_config(u)

source=Source("Wall", source_position=source_position)


#--- Propagator calculation
#
nc=Common.Nc
nd=Common.Nd

print("Solving quark propagator")
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

#--- 2pt correlator

print("2-point correlator:")

corr=Corr2pt_4spinor(GammaMatrixSet("Chiral"))

corr.meson_all(pp, pp)
