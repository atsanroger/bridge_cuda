import math
from pybridge import *

setup(lattice_size=[4,4,4,8], grid_size=[], verbose_level="Detailed")

nvol=Common.Nvol
ndim=Common.Ndim
nc=Common.Nc

u=Field_G(nvol, ndim)

rand=RandomNumbers_Mseries(seed=1234567)
u.set_random(rand)

# check unitarity

av=0.0
for isite in range(nvol):
    for mu in range(ndim):
        w=u.mat(isite,mu)

        # | u u^dag - 1 |
        v=Mat_SU_N(nc)
        v.unit()
        v*=-1.0
        v.multadd_nd(w, w)

        av+=math.sqrt(v.norm2())

print()
print("Random SU({}):".format(nc))
print("  ave |UU^dag - I| = {:23.16e}".format(av/nvol/ndim))
print()
