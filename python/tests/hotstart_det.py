import numpy as np
import scipy.linalg as la
import cmath
from pybridge import *

setup(lattice_size=[4,4,4,8], grid_size=[], verbose_level="Detailed")

nvol=Common.Nvol
ndim=Common.Ndim
nc=Common.Nc

u=Field_G(nvol, ndim)

rand=RandomNumbers_Mseries(seed=1234567)
u.set_random(rand)

uu=np.array(u)

one=np.identity(nc,dtype=complex)


# check unitarity

av=0.0
for isite in range(nvol):
    for mu in range(ndim):

        ww=uu[:,isite,mu].reshape(3,3)

        # | u u^dag - 1 |
        v=la.norm( ww.dot(ww.conj().T) - one )
        
        av+=v

print()
print("Random SU({}): unitarity".format(nc))
print("  ave |UU^dag - I|  = {:23.16e}".format(av/nvol/ndim))
print()


# determinant

re_av=0.0
re_va=0.0
im_av=0.0
im_va=0.0

for isite in range(nvol):
    for mu in range(ndim):

        ww=uu[:,isite,mu].reshape(3,3)

        v=la.det(ww)

        re_av += v.real
        re_va += v.real ** 2
        im_av += v.imag
        im_va += v.imag ** 2
        
re_avg = re_av / nvol / ndim
re_var = re_va / nvol / ndim - re_avg ** 2
im_avg = im_av / nvol / ndim
im_var = im_va / nvol / ndim - im_avg ** 2

print("Random SU({}): determinant".format(nc))
print("  ave Re(det)       = {:23.16e}".format(re_avg))
print("  var Re(det)       = {:23.16e}".format(np.sqrt(re_var)))
print("  ave Im(det)       = {:23.16e}".format(im_avg))
print("  var Im(det)       = {:23.16e}".format(np.sqrt(im_var)))
print()

# determinant using lu decomp

re_av=0.0
re_va=0.0
im_av=0.0
im_va=0.0

for isite in range(nvol):
    for mu in range(ndim):

        ww=uu[:,isite,mu].reshape(3,3)

        P,L,U=la.lu(ww)
        v = np.prod(np.diag(L)) * np.prod(np.diag(U)) * la.det(P)
        
        re_av += v.real
        re_va += v.real ** 2
        im_av += v.imag
        im_va += v.imag ** 2
        
re_avg = re_av / nvol / ndim
re_var = re_va / nvol / ndim - re_avg ** 2
im_avg = im_av / nvol / ndim
im_var = im_va / nvol / ndim - im_avg ** 2

print("Random SU({}): determinant by LU".format(nc))
print("  ave Re(det)       = {:23.16e}".format(re_avg))
print("  var Re(det)       = {:23.16e}".format(np.sqrt(re_var)))
print("  ave Im(det)       = {:23.16e}".format(im_avg))
print("  var Im(det)       = {:23.16e}".format(np.sqrt(im_var)))
print()


# eigenvalues

av=0.0
va=0.0

for isite in range(nvol):
    for mu in range(ndim):

        ww=uu[:,isite,mu].reshape(3,3)

        e,v=la.eig(ww)
        for ic in range(nc):
            x=cmath.phase(e[ic])
            av += x
            va += x ** 2
            
avg=av/nvol/ndim/nc
var=va/nvol/ndim/nc - avg**2

print("Random SU({}): eigenvalues".format(nc))
print("  ave              = {:23.16e}".format(avg))
print("  var              = {:23.16e}".format(np.sqrt(var)))
print("  var(expected)    = {:23.16e}".format(np.pi/np.sqrt(3.0)))
print()
