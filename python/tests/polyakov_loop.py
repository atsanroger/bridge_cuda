import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from pybridge import *

setup(lattice_size=[4,4,4,8], grid_size=[],
      random_number_type="Mseries",
      random_number_seed=1234567)

nvol=Common.Nvol
ndim=Common.Ndim

u=Field_G(nvol, ndim)

io=GaugeConfig("Text")
io.read_file(u, "conf_04040408.txt")

ploop=PolyakovLoop()

result=ploop.measure(u)

print("result=", result)
