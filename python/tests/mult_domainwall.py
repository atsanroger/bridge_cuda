import timeit
from pybridge import *

setup(lattice_size=[4,4,4,8], grid_size=[], verbose_level="Detailed")

nvol=Common.Nvol
ndim=Common.Ndim

nc=Common.Nc
npe=Common.NPE
nthr=Common.Nthread

ns=8  # extent of 5th dimension


u=Field_G(nvol, ndim)

gconf=GaugeConfig("Text")
gconf.read_file(u, "conf_04040408.txt")

nmult=1000

print("--------------------------------")
print("mult performance: nmult={}".format(nmult))
print("--------------------------------")
print()

fopr=Fopr_Domainwall(Fopr_Wilson(gamma_matrix_type="Dirac",
                                 hopping_parameter=1.0/(8.0-2.0*1.6), #dummy
                                 boundary_condition=[1,1,1,-1]),
                     quark_mass=0.1,
                     quark_mass_PauliVillars=1.0,
                     domain_wall_height=1.6,
                     extent_of_5th_dimension=ns,
                     boundary_condition=[1,1,1,-1])

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


# floating-point operations(flo) of mult_Domainwall(Imp version)
#   (D_DW) = (464 * Nc + 240) * Ns * site
#   (D_W) = 464 * Nc * site

nops = 1.0e-9 * (464 * nc + 240) * ns * nvol * npe
gflops = nops * nmult / (elapsed * npe * nthr)

print("elapsed time = {:.6f} sec".format(elapsed))
print("performance  = {:.6f} GFlops".format(gflops))
print()
