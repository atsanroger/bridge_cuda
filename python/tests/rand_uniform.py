from pybridge import *
import math

setup(lattice_size=[4,4,4,8], grid_size=[], verbose_level="Detailed")

nvol=Common.Nvol
ndim=Common.Ndim

nseed=20
seed_base=200000
nsample=10000000

def run_test(type):
    print("Random number type = {}".format(type))
    print("  Nseed      = {}".format(nseed))
    print("  iseed_base = {}".format(seed_base))
    print("  Nrand      = {}".format(nsample))
    print()
    print("Monte Carlo estimate of pi:");
    print("  number of samples = {:10d}".format(nsample))
    print("        seed    estimate of pi")

    t1=0.0
    t2=0.0

    for iseed in range(nseed):
        seed=seed_base+iseed

        rand=RandomNumbers(type, seed=seed)

        npi=0
        for i in range(nsample):
            r1=rand.get()
            r2=rand.get()
            r=r1**2 + r2**2
            if r < 1.0:
                npi += 1

        vpi=4.0 * npi / nsample

        t1+=vpi
        t2+=vpi**2

        print("  {:10d}    {:14.10f}".format(seed, vpi))

    api=t1/nseed
    vpi=t2/nseed - api**2
    dpi=math.sqrt(vpi)
    epi=dpi/math.sqrt(nseed-1)

    print("  true value = {:10.8f}".format(math.pi))
    print("  average    = {:10.8f}".format(api))
    print("  variance   = {:10.8f}".format(vpi))
    print("  deviation  = {:10.8f}".format(dpi))
    print("  error      = {:10.8f}".format(epi))
    print()


run_test("Mseries")
run_test("MT19937")
run_test("SFMT")
