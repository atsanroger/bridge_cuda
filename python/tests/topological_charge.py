from pybridge import *

setup(lattice_size=[4,4,4,8], grid_size=[], verbose_level="Detailed")

nvol=Common.Nvol
ndim=Common.Ndim

u=Field_G(nvol, ndim)

gconf=GaugeConfig("Text")
gconf.read_file(u, "conf_04040408.txt")

smear=Smear("APE",
            Projection("Maximum_SU_N",
                       maximum_number_of_iteration=1000,
                       convergence_criterion=1.0e-14),
            rho_uniform=0.1)

topo=TopologicalCharge(c_plaq=5.0/3.0,
                       c_rect=-1.0/12.0,
                       max_momentum=1,
                       filename_output="stdout")

nsmear=30
nmeas=1

for i in range(nsmear+1):
    print("smear step ", i)

    u2=u.clone()
    if i == 0:
        copy(u2,u)
    else:
        smear.smear(u2,u)

    if i % nmeas == 0:
        r=topo.measure(u2)
        print("i_smear,Q = {} {:20.16e}".format(i,r))

    copy(u,u2)
