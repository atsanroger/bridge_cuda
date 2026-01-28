import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from pybridge import *

setup(lattice_size=[4,4,4,8], grid_size=[],
      random_number_type="Mseries",
      random_number_seed=1234567,
      verbose_level="Detailed")

nvol=Common.Nvol
ndim=Common.Ndim

u=Field_G(nvol, ndim)

gconf=GaugeConfig("Text")
gconf.read_file(u, "conf_04040408.txt")

smear=Smear("APE_spatial",
            Projection("Maximum_SU_N",
                       maximum_number_of_iteration=1000,
                       convergence_criterion=1.0e-14),
            rho=1.0/(1.0+2.3))

wloop=WilsonLoop(max_spatial_loop_size=2,
                 max_temporal_loop_size=4,
                 number_of_loop_type=6,
                 filename_output="stdout",
                 verbose_level="General")

nsmear=30
nmeas=5

for i in range(nsmear+1):
    print("smear step ", i)

    u2=u.clone()
    if i == 0:
        copy(u2,u)
    else:
        smear.smear(u2,u)

    if i % nmeas == 0:
        r=wloop.measure(u2)
        print("result=",r)

    copy(u,u2)
