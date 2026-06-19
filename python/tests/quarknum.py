from pybridge import *

setup(lattice_size=[4,4,4,8], grid_size=[], verbose_level="Detailed")

nvol=Common.Nvol
ndim=Common.Ndim

u=Field_G(nvol, ndim)

gconf=GaugeConfig("Text")
gconf.read_file(u, "conf_04040408.txt")

smear=Director_Smear(Smear("APE",
                           Projection("Stout_SU3"),
                           rho_uniform=0.1),
                     number_of_smearing=1)

fopr=Fopr_Smeared(Fopr_Clover_Isochemical(gamma_matrix_type="Dirac",
                                          hopping_parameter=0.12,
                                          clover_coefficient=1.0,
                                          isospin_chemical_potential=0.3,
                                          boundary_condition=[1,1,1,-1]),
                  smear)
fopr.set_config(u)

solver=Solver(type="BiCGStab_L_Cmplx",
              op=fopr,
              maximum_number_of_iteration=100,
              maximum_number_of_restart=40,
              convergence_criterion_squared=1.0e-28,
              use_initial_guess="False",
              Omega_tolerance=0.60,
              number_of_orthonormal_vectors=2)
fprop=Fprop_Standard_lex(solver)

rand=RandomNumbers("Mseries", seed=200000)
noise=NoiseVector_Z2(rand)

qsuscept=QuarkNumberSusceptibility_Wilson(fopr, fprop, noise,
                                          number_of_noises=10)

result=qsuscept.measure()

print("result=", result)
