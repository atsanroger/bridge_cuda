import timeit
from pybridge import *

setup(lattice_size=[4,4,4,8], grid_size=[], verbose_level="Detailed")

nvol=Common.Nvol
ndim=Common.Ndim

u=Field_G(nvol, ndim)

gconf=GaugeConfig("Text")
gconf.read_file(u, "conf_04040408.txt")

nmult=1000

print("--------------------------------")
print("mult performance: nmult={}".format(nmult))
print("--------------------------------")
print()


def run_test(name, fopr):
    print("---- {} ----".format(name))

    fopr.set_mode("D")
    fopr.set_config(u)

    b=Field_F(nvol)

    nin=b.nin
    nex=b.nex
    
    be=Field(nin, nvol//2, nex)
    ye=Field(nin, nvol//2, nex)
    be.fill(1.0)

    elapsed=timeit.timeit(
        lambda:
          fopr.mult(ye, be),
        number=nmult)

    result=ye.norm()

    gflops=fopr.flop_count() * nmult / elapsed

    print("elapsed time = {:.6f} sec".format(elapsed))
    print("performance  = {:.6f} GFlops".format(gflops))
    print("-" * (len(name) + 10))
    print()


run_test("Fopr_Wilson_eo",
         Fopr("Wilson_eo",
              gamma_matrix_type="Dirac",
              hopping_parameter=0.12,
              boundary_condition=[-1,-1,-1,-1])
)

run_test("Fopr_Clover_eo",
         Fopr("Clover_eo",
              gamma_matrix_type="Dirac",
              hopping_parameter=0.12,
              clover_coefficient=1.0,
              boundary_condition=[1,1,1,-1])
)
