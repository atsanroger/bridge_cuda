from pybridge import *
import math

setup(lattice_size=[4,4,4,8], grid_size=[], verbose_level="Detailed")

nvol=Common.Nvol
ndim=Common.Ndim

nin=30
nex=33
size=nin*nvol*nex


def run_test(rand, type):
    print("  noise_type = {}".format(type))
    
    f=Field(nin,nvol,nex)

    rand.reset(0)
    rand.lex_global(type, f)

    av=0.0
    vr=0.0

    for i in range(size):
        x=f.cmp(i)
        av += x
        vr += x * x

    av=av/nvol/nin/nex
    vr=vr/nvol/nin/nex - av**2
    vr=math.sqrt(vr)

    print("Distribution (Field):")
    print("  number of samples = {:10d}".format(size))
    print("  average           = {:10.8f}".format(av))
    print("  variance          = {:10.8f}".format(vr))
    print("  variance(expect)  = {:10.8f}".format(math.sqrt(0.5)))
    print()


print("RandomNumbers: Mseries")
rand=RandomNumbers_Mseries(seed=0)

if rand:
    run_test(rand, "Gaussian")
    run_test(rand, "U1")
    run_test(rand, "Z2")


print("RandomNumbers: MT19937")
rand=RandomNumbers_MT19937(seed=0)

if rand:
    run_test(rand, "Gaussian")
    run_test(rand, "U1")
    run_test(rand, "Z2")


print("RandomNumbers: SFMT")
rand=RandomNumbers_SFMT(seed=0)

if rand:
    run_test(rand, "Gaussian")
    run_test(rand, "U1")
    run_test(rand, "Z2")
