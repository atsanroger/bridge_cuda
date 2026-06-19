from pybridge import *

setup(lattice_size=[4,4,4,8], grid_size=[], verbose_level="Detailed")

nvol=Common.Nvol
ndim=Common.Ndim

u=Field_G(nvol, ndim)

gconf=GaugeConfig("Text")
gconf.read_file(u, "conf_04040408.txt")

#gfix=GaugeFixing("None")
gfix=GaugeFixing("Coulomb",
                 maximum_number_of_iteration=5000,
                 number_of_naive_iteration=50,
                 interval_of_measurement=10,
                 iteration_to_reset=1000,
                 convergence_criterion_squared=1.0e-24,
                 overrelaxation_parameter=1.6,
                 verbose_level="Detailed"
)

ufix=u.clone()
gfix.fix(ufix,u)
copy(u,ufix)


source=Source("Local",
              source_position=[0,0,0,0],
)

fopr=Fopr("Wilson",
          gamma_matrix_type="Dirac",
          hopping_parameter=0.12,
          boundary_condition=[-1,-1,-1,-1])

fopr.set_config(u)


def test(solver_type, solver):
    print("--- {} ---".format(solver_type))
    
    fprop=Fprop_Standard_lex(solver)

    nc=Common.Nc
    nd=Common.Nd

    for ispin in range(nd):
        for icolor in range(nc):
            b=Field_F(nvol)
            source.set(b, icolor=icolor, ispin=ispin)

            x=Field_F(nvol)
            nconv,diff=fprop.invert_D(x, b)

            y=Field_F(nvol)
            copy(y,b)
        
            fopr.set_mode("D")
            fopr.mult(y,x)
            axpy(y, -1.0, b)
            diff2=y.norm2()
        
            print("   {:2d}   {:2d}   {:6d}   {:12.4e}   {:12.4e}".format(icolor,ispin,nconv,diff,diff2))


test("BiCGStab_Cmplx",
    Solver(type="BiCGStab_Cmplx",
           op=fopr,
           maximum_number_of_iteration=100,
           convergence_criterion_squared=1.0e-28,
           maximum_number_of_restart=40,
           Omega_tolerance=0.7,
           use_initial_guess="false"
           ))

test("BiCGStab_L_Cmplx",
    Solver(type="BiCGStab_L_Cmplx",
           op=fopr,
           maximum_number_of_iteration=100,
           convergence_criterion_squared=1.0e-28,
           maximum_number_of_restart=40,
           Omega_tolerance=0.7,
           number_of_orthonormal_vectors=2,
           use_initial_guess="false"
           ))

test("BiCGStab_DS_L_Cmplx",
    Solver(type="BiCGStab_DS_L_Cmplx",
           op=fopr,
           maximum_number_of_iteration=100,
           convergence_criterion_squared=1.0e-28,
           maximum_number_of_restart=40,
           Omega_tolerance=0.7,
           number_of_orthonormal_vectors=2,
           tolerance_for_DynamicSelection_of_L=0.1,
           use_initial_guess="false"
           ))

test("BiCGStab_IDS_L_Cmplx",
    Solver(type="BiCGStab_IDS_L_Cmplx",
           op=fopr,
           maximum_number_of_iteration=100,
           convergence_criterion_squared=1.0e-28,
           maximum_number_of_restart=40,
           Omega_tolerance=0.7,
           number_of_orthonormal_vectors=2,
           tolerance_for_DynamicSelection_of_L=0.1,
           use_initial_guess="false"
           ))

test("CGNE",
    Solver(type="CGNE",
           op=fopr,
           maximum_number_of_iteration=100,
           convergence_criterion_squared=1.0e-28,
           maximum_number_of_restart=40,
           use_initial_guess="false"
           ))

test("CGNR",
    Solver(type="CGNR",
           op=fopr,
           maximum_number_of_iteration=100,
           convergence_criterion_squared=1.0e-28,
           maximum_number_of_restart=40,
           use_initial_guess="false"
           ))

test("GMRES_m_Cmplx",
    Solver(type="GMRES_m_Cmplx",
           op=fopr,
           maximum_number_of_iteration=100,
           convergence_criterion_squared=1.0e-28,
           maximum_number_of_restart=40,
           number_of_orthonormal_vectors=13,
           use_initial_guess="false"
           ))

