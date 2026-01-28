from pybridge import *
import numpy as np

setup(lattice_size=[4,4,4,8], grid_size=[], verbose_level="Detailed")

nvol=Common.Nvol
ndim=Common.Ndim

#--- local volume and data layout
layout=(Common.Nc**2, Common.Nx, Common.Ny, Common.Nz, Common.Nt, ndim)

u=Field_G(nvol,ndim)

gconf=GaugeConfig("Text")
gconf.read_file(u, "conf_04040408.txt")

shift=ShiftField_lex()

#--- original field
print("original field:")

Staple_lex().plaquette(u)

# export to numpy array
uu=np.array(u).reshape(layout, order='F')

#--- shift backward along t-direction
print("shift backward along t-direction:")
u2=u.clone()
shift.backward(u2, u, 3)

Staple_lex().plaquette(u2)

uu2=np.array(u2).reshape(layout, order='F')
print("before backward=", uu[0, 0,0,0,:, 0])
print("after backward= ", uu2[0, 0,0,0,:, 0])

if Common.NPEt == 1:
    print("check", np.array_equal(np.roll(uu, -1, axis=4), uu2))

print()

#--- shift forward along x-direction
print("shift forward along x-direction:")
u3=u.clone()
shift.forward(u3, u, 0)

Staple_lex().plaquette(u3)

uu3=np.array(u3).reshape(layout, order='F')
print("before forward= ", uu[0, :,0,0,0, 0])
print("after forward=  ", uu3[0, :,0,0,0, 0])

if Common.NPEx == 1:
    print("check", np.array_equal(np.roll(uu, 1, axis=1), uu3))

print()

