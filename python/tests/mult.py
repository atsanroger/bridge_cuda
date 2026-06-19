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
    print("-" * (len(name) + 10))
    print()


run_test("Fopr_Wilson",
         Fopr("Wilson",
              gamma_matrix_type="Dirac",
              hopping_parameter=0.12,
              boundary_condition=[-1,-1,-1,-1])
)

run_test("Fopr_WilsonGeneral",
         Fopr("WilsonGeneral",
              gamma_matrix_type="Dirac",
              hopping_parameter_spatial=0.10959947,
              hopping_parameter_temporal=0.10959947,
              dispersion_parameter_spatial=1.1450511,
              Wilson_parameter_spatial=1.1881607,
              boundary_condition=[-1,-1,-1,-1])
)

run_test("Fopr_Wilson_Isochemical",
         Fopr("Wilson_Isochemical",
              gamma_matrix_type="Dirac",
              hopping_parameter=0.12,
              isospin_chemical_potential=0.3,
              boundary_condition=[-1,-1,-1,-1])
)

run_test("Fopr_Clover",
         Fopr("Clover",
              gamma_matrix_type="Dirac",
              hopping_parameter=0.12,
              clover_coefficient=1.0,
              boundary_condition=[1,1,1,-1])
)

run_test("Fopr_CloverGeneral",
         Fopr("CloverGeneral",
              gamma_matrix_type="Dirac",
              hopping_parameter_spatial=0.10959947,
              hopping_parameter_temporal=0.10959947,
              dispersion_parameter_spatial=1.1450511,
              Wilson_parameter_spatial=1.1881607,
              clover_coefficient_spatial=1.9849139,
              clover_coefficient_temporal=1.7819512,
              boundary_condition=[1,1,1,-1])
)

run_test("Fopr_Clover_Isochemical",
         Fopr("Clover_Isochemical",
              gamma_matrix_type="Dirac",
              hopping_parameter=0.12,
              clover_coefficient=1.0,
              isospin_chemical_potential=0.3,
              boundary_condition=[-1,-1,-1,-1])
)

