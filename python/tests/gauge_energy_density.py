from pybridge import *

setup(lattice_size=[4,4,4,8], grid_size=[], verbose_level="Detailed")

nvol=Common.Nvol
ndim=Common.Ndim

u=Field_G(nvol,ndim)

gconf=GaugeConfig("Text")
gconf.read_file(u, "conf_04040408.txt")


energy_density=EnergyDensity(c_plaq=5.0/3.0,
                             c_rect=-1.0/12.0)

energy_density.E_clover(u)
