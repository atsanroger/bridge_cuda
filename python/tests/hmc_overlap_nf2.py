from pybridge import *

setup(lattice_size=[4,4,4,8], grid_size=[], verbose_level="General")

nvol=Common.Nvol
ndim=Common.Ndim

#-- gauge part

action_g=Action_G_Rectangle(beta=2.3,
                            c_plaq=3.648,
                            c_rect=-0.331)

#-- Nf=2 part

fopr_params=dict(gamma_matrix_type="Dirac",
                 hopping_parameter=0.5/(4.0-1.6),
                 quark_mass=0.2,
                 domain_wall_height=1.6,
                 number_of_poles=16,
                 lower_bound=0.01,
                 upper_bound=8.0,
                 maximum_number_of_iteration=1000,
                 convergence_criterion_squared=1.0e-20,
                 boundary_condition=[1,1,1,1])

fopr_w=Fopr_Wilson(**fopr_params)

fopr=Fopr_Overlap(fopr_w, **fopr_params)
force_fopr=Force_F_Overlap_Nf2(**fopr_params)

solver_md=Solver_CG(fopr,
                    maximum_number_of_iteration=100,
                    maximum_number_of_restart=40,
                    convergence_criterion_squared=1.0e-24,
                    use_initial_guess="false")
fprop_md=Fprop_Standard_lex(solver_md)

solver_h=Solver_CG(fopr,
                   maximum_number_of_iteration=100,
                   maximum_number_of_restart=40,
                   convergence_criterion_squared=1.0e-24,
                   use_initial_guess="false")
fprop_h=Fprop_Standard_lex(solver_h)


action_f=Action_F_Overlap_Nf2(fopr, force_fopr, fprop_md, fprop_h)


#-- integrator

actions=[action_f, action_g]

directors=[]

#builder=Builder_Integrator(actions, directors,
#                           integrator="Leapfrog",
#                           #number_of_levels=1,
#                           number_of_steps=[5],
#                           order_of_exp_iP=8,
#                           lambda_Omelyan=0.1931833275)
#integ=builder.build()


#-- hmc
traj_number=200000

rand=RandomNumbers_Mseries(traj_number)

#hmc=HMC_General(actions, directors, integ, rand,
#                trajectory_length=0.1,
#                Metropolis_test="true")
hmc=HMC_Leapfrog(actions, directors, rand,
                 step_size=0.02,
                 number_of_steps=5,
                 order_of_exp_iP=8,
                 Metropolis_test="true")

#-- exec

u=Field_G(nvol,ndim)

gconf=GaugeConfig("Text")
gconf.read_file(u, "conf_04040408.txt")

gconfout=GaugeConfig("None") # for output


ntraj=1
ntraj_interval=100

for traj in range(ntraj):
    print("")
    print("traj = {}".format(traj))

    traj_id = traj_number + traj + 1

    r=hmc.update(u)

    if traj_id % ntraj_interval == 0:
        print("---- write config ----")
        gconfout.write_file(u, "conf_output.dat-{:06}".format(traj_id))

gconfout.write_file(u, "conf_output.dat")

print("---- done ----")
