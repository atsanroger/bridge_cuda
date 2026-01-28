from pybridge import *

setup(lattice_size=[4,4,4,8], grid_size=[],
      random_number_type="Mseries",
      random_number_seed=1234567,
      verbose_level="Detailed")

nvol=Common.Nvol
ndim=Common.Ndim

u=Field_G(nvol,ndim)

gconf=GaugeConfig("Text")
gconf.read_file(u, "conf_04040408.txt")

smear=Director_Smear(Smear("APE",
                           Projection("Stout_SU3"),
                           rho_uniform=0.1),
                     number_of_smearing=0)

fopr_core=Fopr_Smeared(Fopr("Clover",
                            gamma_matrix_type="Dirac",
                            hopping_parameter=0.12,
                            clover_coefficient=1.0,
                            boundary_condition=[1,1,1,-1]),
                       smear)
fopr=Fopr_Chebyshev(fopr_core,
                    degree_of_polynomial=40,
                    threshold_value=0.16,
                    upper_bound=2.50)
fopr_core.set_config(u)
fopr_core.set_mode("DdagD")
#XXX workaround for unimplemented set_mode in Fopr_Chebyshev.

vthr=fopr.mult(0.15**2)
print("Vthreshold_ch = {:12.6f}".format(vthr))

# N.B. Low and High must be inverted, because of Chebyshev expansion.
eigen=Eigensolver_IRLanczos(fopr,
                            eigensolver_mode="abs_descending",
                            number_of_wanted_eigenvectors=20,
                            number_of_working_eigenvectors=50,
                            convergence_criterion_squared=1.0e-26,
                            maximum_number_of_iteration=500,
                            threshold_value=vthr,
                            verbose_level="Detailed")

b=Field_F(nvol)

vk,tda,nsbt,nconv=eigen.solve(b,work_size=70)

for i in range(nsbt+1):
    v=Field_F(nvol)

    fopr_core.set_mode("H")
    fopr_core.mult(v,vk[i])

    vnum=dot(vk[i],v)
    vden=dot(vk[i],vk[i])
    veig=vnum/vden
    
    axpy(v, -veig, vk[i])
    vv=v.norm2()

    # print("Eigenvalues: {} {:20.14f} {:20.15e}".format(i,tda[i],vv))
    print("Eigenvalues: {:4d} {:20.14f}  {:10.4e}  {:10.4e}"
          .format(i,veig,vv, vden-1.0))
