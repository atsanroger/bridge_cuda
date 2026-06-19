from pybridge import *

setup(lattice_size=[4,4,4,8], grid_size=[], verbose_level="Detailed")

nvol=Common.Nvol
ndim=Common.Ndim

if Common.NPE > 1:
    print("spectrum_CRS: parallelization not supported. skip.")
    exit()

#--- Parameters

matrix_file="matrix4x8_clover.crs"
source_file="source4x8_clover.crs"
solution_file="solution4x8_clover.crs"

#--- Load gauge config

u=Field_G(nvol, ndim)

gconf=GaugeConfig("Text")
gconf.read_file(u, "conf_04040408.txt")

#--- Setup objects

fopr_c=Fopr_Clover(gamma_matrix_type="Dirac",
                   hopping_parameter=0.12,
                   clover_coefficient=1.0,
                   boundary_condition=[-1,-1,-1,-1])
fopr_c.set_config(u)
fopr_c.set_mode("D")


fopr_crs=Fopr_CRS(fopr_c)
fopr_crs.write_matrix(matrix_file)


solver=Solver_CG(fopr_c,
                 maximum_number_of_iteration=100,
                 maximum_number_of_restart=40,
                 convergence_criterion_squared=1.0e-28,
                 use_initial_guess="false")

source=Source_Local(source_position=[0,0,0,0])

#--- tools

def write_text(v, filename):
    print("Writing Field to {}".format(filename))

    with open(filename, mode="w") as f:
        print(v.nin, file=f)
        print(v.nvol, file=f)
        print(v.nex, file=f)

        for j in range(v.size):
            print(format(v.cmp(j), ".14e"), file=f)

    print("Writing Field finished.")

def read_text(v, filename):
    print("Reading Field from {}".format(filename))

    try:
        with open(filename, mode="r") as f:
            nin=int(f.readline())
            nvol=int(f.readline())
            nex=int(f.readline())

            v.reset(nin, nvol, nex)

            for j in range(v.size):
                v.set(j, float(f.readline()))

    except FileExistsError:
        print("Error: Failed to open the text file. {}".format(filename))

    print("Reading Field finished.")

#--- Propagator calculation

nc=Common.Nc
nd=Common.Nd

pp=[0] * (nc * nd)

for ispin in range(nd):
    for icolor in range(nc):
        icd=icolor + nc * ispin
        
        b=Field_F(nvol)
        source.set(b, icolor=icolor, ispin=ispin)

        write_text(b, source_file)
        
        b2=b.clone()
        fopr_c.set_mode("D")
        fopr_c.mult_dag(b2, b)
        
        xq=b.clone()

        fopr_c.set_mode("DdagD")
        nconv, diff = solver.solve(xq, b2)

        print("  ispin = {:2d}  icolor = {:2d}  Nconv = {:4d}  diff = {:12.6e}"
              .format(ispin, icolor, nconv, diff))

        write_text(xq, solution_file)

        y=b.clone()
        copy(y, b)

        fopr_c.set_mode("D")
        fopr_c.mult(y, xq)
        axpy(y, -1.0, b)
        scal(y, -1.0)

        diff2=y.norm2()

        print("    standard norm2 = {:.8e}".format(diff2))

        pp[icd] = xq

print(pp)


def crs_solver(solution, matrix, source):
    b=Field()
    read_text(b, source)

    xq=Field()
    read_text(xq, solution)

    fopr_crs=Fopr_CRS(matrix)

    solver=Solver_CG(fopr_crs,
                     maximum_number_of_iteration=100,
                     maximum_number_of_restart=40,
                     convergence_criterion_squared=1.0e-28,
                     use_initial_guess="false")

    # setup CGNE source
    b2=b.clone()
    fopr_crs.set_mode("D")
    fopr_crs.mult_dag(b2, b)

    # CGNE solver
    x=b.clone()
    fopr_crs.set_mode("DdagD")
    nconv,diff = solver.solve(x, b2)
    print("solver(CG):  Nconv = {:4d}  diff = {:12.6e}".format(nconv, diff))

    # check
    axpy(x, -1.0, xq)
    xx=x.norm2()

    print("standard norm2 = {:.8e}".format(xx))

    return xx


result=crs_solver(solution_file, matrix_file, source_file)

print("result = ", result)


