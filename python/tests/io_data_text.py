from pybridge import *
import numpy as np
import os

setup(lattice_size=[4,4,4,8], grid_size=[])

nvol=Common.Nvol
ndim=Common.Ndim

nx=Common.Nx
ny=Common.Ny
nz=Common.Nz
nt=Common.Nt
nc=Common.Nc

u=Field_G(nvol, ndim)

gconf=GaugeConfig("Text")
gconf.read_file(u, "conf_04040408.txt")

uu=np.array(u).reshape(nc,nc,nx,ny,nz,nt,ndim,order='F')


v=abs(uu[:,:,:,:,:,0,0]).flatten()
#print(v)

datalen=len(v)
print("data length =", datalen)

if Communicator.self() == 0:
    if os.path.isfile("output.dat"):
        print("output file exists. rename.")
        os.rename("output.dat", "output.dat.bak")

dataio=DataIO("Text")
dataio.write_file(v.tolist(), "output.dat", append=False)

v2=np.array(dataio.read_file(file="output.dat", length=datalen))
#print(v2)


if Communicator.self() == 0:
    #tol=1.0e-15
    tol=1.0e-14

    ok=abs((v-v2)/v) < tol

    if np.all(ok) == True:
        print("ok")
    else:
        idx=np.where(ok == False)
        print("failed. err=", np.count_nonzero(idx))
        for x in list(zip(*idx,v[idx],v2[idx]))[0:32]:
            print(*x)
