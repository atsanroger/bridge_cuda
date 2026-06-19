from pybridge import *
import numpy as np
import os

setup(lattice_size=[4,4,4,8], grid_size=[])

nvol=Common.Nvol
ndim=Common.Ndim

u=Field_G(nvol, ndim)

gconf=GaugeConfig("Text")
gconf.read_file(u, "conf_04040408.txt")


def test(conf_type):
    print("--- {} ---".format(conf_type))

    io=GaugeConfig(conf_type)

    filename="output_" + conf_type + ".dat"
    print(filename)
    
    if Communicator.self() == 0:
        if os.path.isfile(filename):
            print("output file exists. rename.")
            os.rename(filename, filename + ".bak")

    io.write_file(u, filename)

    uu=u.clone()
    io.read_file(uu, filename)

    v=np.array(u)
    vv=np.array(uu)

    if np.array_equal(v, vv):
        print("ok")
    else:
        print("failed")

    return

def test_random():
    conf_type="Random"
    print("--- {} ---".format(conf_type))

    nc=Common.Nc
    zero=np.zeros((nc,nc),dtype=complex)
    one=np.identity(nc,dtype=complex)

    io=GaugeConfig(conf_type)

    filename="output_" + conf_type + ".dat"
    print(filename)
    
    uu=u.clone()
    io.read_file(uu, filename)

    vv=np.array(uu).reshape(nc,nc,nvol,ndim,order='F')

    err=0
    for idim in range(ndim):
        for isite in range(nvol):
            w=vv[:,:,isite,idim]
            ww=np.conjugate(w.T)
            x=np.matmul(w,ww)

            if np.allclose(x,one,atol=1.0e-14) == False:
                err += 1
                print(isite, idim, np.linalg.norm(x))

    if err == 0:
        print("ok")
    else:
        print("failed")

    return

def test_unit():
    conf_type="Unit"
    print("--- {} ---".format(conf_type))

    nc=Common.Nc
    zero=np.zeros((nc,nc),dtype=complex)
    one=np.identity(nc,dtype=complex)

    io=GaugeConfig(conf_type)

    filename="output_" + conf_type + ".dat"
    print(filename)
    
    uu=u.clone()
    io.read_file(uu, filename)

    vv=np.array(uu).reshape(nc,nc,nvol,ndim,order='F')

    err=0
    for idim in range(ndim):
        for isite in range(nvol):
            w=vv[:,:,isite,idim]
            if np.array_equal(w,one) == False:
                err+=1
                print("error in {},{}".format(isite,idim))

    if err == 0:
        print("ok")
    else:
        print("failed")

    return

def test_none():
    conf_type="None"
    print("--- {} ---".format(conf_type))

    nc=Common.Nc
    zero=np.zeros((nc,nc),dtype=complex)
    one=np.identity(nc,dtype=complex)

    io=GaugeConfig(conf_type)

    filename="output_" + conf_type + ".dat"
    print(filename)
    
    uu=u.clone()
    io.read_file(uu, filename)

    vv=np.array(uu).reshape(nc,nc,nvol,ndim,order='F')

    err=0
    for idim in range(ndim):
        for isite in range(nvol):
            w=vv[:,:,isite,idim]
            if np.array_equal(w,zero) == False:
                err+=1
                print("error in {},{}".format(isite,idim))

    if err == 0:
        print("ok")
    else:
        print("failed")

    return


test("Text")
test("Binary")
test("Fortran_ILDG")
test("Fortran_JLQCD")

test("ILDG")

if Communicator.mpi_enabled() == True:
    test("Binary_Parallel")
    test("Binary_Distributed")
    test("ILDG_Parallel")

test_random()
test_unit()
test_none()
