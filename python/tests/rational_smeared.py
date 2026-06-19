import math
from pybridge import *

setup(lattice_size=[4,4,4,8], grid_size=[], verbose_level="Detailed")

nvol=Common.Nvol
ndim=Common.Ndim

u=Field_G(nvol, ndim)

#--- Load gauge config

gconf=GaugeConfig("Text")
gconf.read_file(u, "conf_04040408.txt")

#--- Smearing

smear=Director_Smear(Smear("APE",
                           Projection("Stout_SU3",
                                      maximum_number_of_iteration=1000),
                           rho_uniform=0.1),
                     number_of_smearing=2)

#--- Source

source=Source("Local", source_position=[0,0,0,0])

#--- Fopr

rational_params=dict(number_of_poles=16,
                     exponent_numerator=-1,
                     exponent_denominator=2,
                     lower_bound=1.0e-2,
                     upper_bound=math.sqrt(10.0),
                     maximum_number_of_iteration=1000,
                     convergence_criterion_squared=1.0e-24)

fopr_c=Fopr("Clover",
            gamma_matrix_type="Dirac",
            hopping_parameter=0.12,
            clover_coefficient=1.0,
            boundary_condition=[1,1,1,-1])

fopr_s=Fopr_Smeared(fopr_c, smear)

fopr_smeared_rational=Fopr_Smeared(
                          Fopr_Rational(
                              fopr_c,
                              **rational_params
                          ), smear)

fopr_rational_smeared=Fopr_Rational(
                          Fopr_Smeared(
                              fopr_c, smear
                          ), **rational_params)

fopr_s.set_config(u)
fopr_smeared_rational.set_config(u)
fopr_rational_smeared.set_config(u)


print("Consistency check of smeared-rational and rational-smeared.")
print()

nc=Common.Nc
nd=Common.Nd

print("Smeared-rational:")

# for ispin in range(nd):
#     for icolor in range(nc):
for ispin in range(1):
    for icolor in range(1):

        b=Field_F(nvol)
        source.set(b, icolor=icolor, ispin=ispin)
        
        v=b.clone()
        fopr_smeared_rational.mult(v, b)

        w=v.clone()
        fopr_s.set_mode("DdagD")
        fopr_s.mult(w, v)

        fopr_smeared_rational.mult(v, w)

        axpy(v, -1.0, b)

        vv=v.norm2()
        print("  standard norm2 = {:.8e}".format(vv))


print("Rational-smeared:")

# for ispin in range(nd):
#     for icolor in range(nc):
for ispin in range(1):
    for icolor in range(1):

        b=Field_F(nvol)
        source.set(b, icolor=icolor, ispin=ispin)
        
        v=b.clone()
        fopr_rational_smeared.mult(v, b)

        w=v.clone()
        fopr_s.set_mode("DdagD")
        fopr_s.mult(w, v)

        fopr_rational_smeared.mult(v, w)

        axpy(v, -1.0, b)

        vv=v.norm2()
        print("  standard norm2 = {:.8e}".format(vv))


