from pybridge import *

setup(lattice_size=[4,4,4,8], grid_size=[],
      random_number_type="Mseries",
      random_number_seed=1234567,
      verbose_level="Detailed")

nvol=Common.Nvol

def run_test(fft, src, ref):

    res1=src.clone()

    fft.fft(res1, src, FFT_Forward)

    res2=src.clone()

    fft.fft(res2, res1, FFT_Backward)

    print("source norm2  = {:24.20e}".format(src.norm2()))
    print("result norm2  = {:24.20e}".format(res1.norm2()))
    print("reverse norm2 = {:24.20e}".format(res2.norm2()))

    # check forward transformation
    axpy(res1, -1.0, ref)

    diff1 = res1.norm()

    # check forward-backward reverse
    axpy(res2, -1.0, src)

    diff2 = res2.norm()

    print("diff (result - ref)  = {:24.20e}".format(diff1))
    print("diff (reverse - src) = {:24.20e}".format(diff2))
    print()

nex=1

src=Field_F(nvol, nex)
ref=src.clone()

NPEx = Common.NPEx
NPEy = Common.NPEy
NPEz = Common.NPEz
NPEt = Common.NPEt

# generate random source and reference

rand=RandomNumbers_MT19937(seed=0xbeef)
rand.gauss_lex_global(src)

fft=FFT("auto")
fft.fft(ref, src, FFT_Forward)


# FFT_xyz_1dim
#
if NPEt == 1 and (NPEx * NPEy == 1 or NPEx * NPEz == 1 or NPEy * NPEz == 1):
    print("run FFT_xyz_1dim")
    run_test(
        FFT("FFT_xyz_1dim", FFT_direction="Forward"),
        src, ref)
else:
    print("skip FFT_xyz_1dim")

print()

# FFT_xyz_3dim
#
if NPEt == 1 and NPEx == 1 and NPEy == 1:
    print("run FFT_xyz_3dim")
    run_test(
        FFT("FFT_xyz_3dim", FFT_direction="Forward"),
        src, ref)
else:
    print("skip FFT_xyz_3dim")

print()

# FFT_3d_local
#
if NPEx * NPEy * NPEz == 1:
    print("run FFT_3d_local")
    run_test(
        FFT("FFT_3d_local"),
        src, ref)
else:
    print("skip FFT_3d_local")

print()

# FFT_3d_parallel1d
#
if Communicator.mpi_enabled() == True:
    print("run FFT_3d_parallel1d")
    if NPEx * NPEy == 1:
        run_test(
            FFT("FFT_3d_parallel1d"),
            src, ref)
    else:
        print("skip FFT_3d_parallel1d")

print()

# FFT_3d_parallel3d
#
if Communicator.mpi_enabled() == True:
    print("run FFT_3d_parallel3d")
    run_test(
        FFT("FFT_3d_parallel3d"),
        src, ref)

print()    

# FFT_3d auto
#
print("run auto")
run_test(
    FFT("auto"),
    src, ref)

