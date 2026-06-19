from pybridge import *

setup(lattice_size=[4,4,4,8], grid_size=[], verbose_level="Detailed")

nvol=Common.Nvol
ndim=Common.Ndim

u=Field_G(nvol,ndim)

gconf=GaugeConfig("Text")
gconf.read_file(u, "conf_04040408.txt")

smear=Director_Smear(Smear("APE",
                           Projection("Stout_SU3"),
                           rho_uniform=0.1),
                     number_of_smearing=0)

fopr=Fopr_Smeared(Fopr("Clover",
                       gamma_matrix_type="Dirac",
                       hopping_parameter=0.12,
                       clover_coefficient=1.0,
                       boundary_condition=[1,1,1,-1]),
                  smear)
fopr.set_config(u)
fopr.set_mode("H")

eigen=Eigensolver_IRLanczos(fopr,
                            eigensolver_mode="abs_ascending",
                            number_of_wanted_eigenvectors=20,
                            number_of_working_eigenvectors=50,
                            convergence_criterion_squared=1.0e-24,
                            maximum_number_of_iteration=500,
                            threshold_value=0.15)

b=Field_F(nvol)

vk,tda,nsbt,nconv=eigen.solve(b,work_size=70)

for i in range(nsbt+1):
    v=Field_F(nvol)

    fopr.mult(v,vk[i])
    axpy(v, -tda[i], vk[i])
    vv=v.norm2()

    print("Eigenvalues: {:4d} {:20.14f} {:20.15e} ".format(i,tda[i],vv))
