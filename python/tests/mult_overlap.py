import timeit
from pybridge import *

setup(lattice_size=[4,4,4,8], grid_size=[], verbose_level="General")

nvol=Common.Nvol
ndim=Common.Ndim

u=Field_G(nvol, ndim)

gconf=GaugeConfig("Text")
gconf.read_file(u, "conf_04040408.txt")

nmult=100

print("--------------------------------")
print("mult performance: nmult={}".format(nmult))
print("--------------------------------")
print()

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

fopr.set_mode("D")
fopr.set_config(u)

nin_dw=fopr.field_nin()
nvol_dw=fopr.field_nvol()
nex_dw=fopr.field_nex()

b=Field(nin_dw, nvol_dw, nex_dw)
b.fill(1.0)

y=b.clone()

elapsed=timeit.timeit(
    lambda:
    fopr.mult(y, b),
    number=nmult)

result=y.norm()

gflops=fopr.flop_count() * nmult / elapsed

print("elapsed time = {:.6f} sec".format(elapsed))
print("performance  = {:.6f} GFlops".format(gflops))
print()
