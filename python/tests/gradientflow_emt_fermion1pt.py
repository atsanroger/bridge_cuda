from pybridge import *

setup(lattice_size=[4,4,4,8], grid_size=[], verbose_level="Detailed")

nvol=Common.Nvol
ndim=Common.Ndim

u=Field_G(nvol,ndim)

gconf=GaugeConfig("Text")
gconf.read_file(u, "conf_04040408.txt")


action=Action_G_Rectangle(c_plaq=1.0, c_rect=0.0, beta=float(Common.Nc))


params_measure=dict(filename_output="stdout",
                    number_of_noises=2,
                    initial_tau=1,
                    number_of_measurement_times=2,
                    measurement_interval=5,
                    gauge_store_interval=1,
                    max_momentum=1)

params_gflow=dict(order_of_RungeKutta=3,
                  step_size=0.02,
                  order_of_approx_for_exp_iP=8,
                  adaptive="false",
                  tolerance=1.0e-4,
                  safety_factor=0.9)

params_fflow=dict(order_of_RungeKutta=3,
                  step_size=0.02,
                  order_of_approx_for_exp_iP=8,
                  boundary_condition=[1,1,1,-1])

params_source=dict(noise_type="Z2",
                   source_position=[0,0,0,0],
                   source_momentum=[0,0,0,0])

params_fopr=[
    dict(type="Clover_eo",
         gamma_matrix_type="Dirac",
         hopping_parameter=0.12,
         clover_coefficient=1.0,
         boundary_condition=[1,1,1,-1]),
    dict(type="Clover_eo",
         gamma_matrix_type="Dirac",
         hopping_parameter=0.115,
         clover_coefficient=1.0,
         boundary_condition=[1,1,1,-1])
]

params_solver=dict(type="BiCGStab_L_Cmplx",
                   maximum_number_of_iteration=5000,
                   maximum_number_of_restart=100,
                   convergence_criterion_squared=1.0e-28,
                   Omega_tolerance=0.60,
                   number_of_orthonormal_vectors=2,
                   use_initial_guess="false")

fprops=[
    Fprop_Standard_eo(
        Solver(
            op=Fopr(**params), **params_solver
        )
    )
    for params in params_fopr
]
print(fprops)

for p in fprops:
    p.set_config(u)

meas=FermionFlow_1pt_Function(fprops, action, params_fopr,
                              params_measure,
                              params_gflow,
                              params_fflow,
                              params_source)

meas.measure_disconnected(u)

# TODO: save random state

