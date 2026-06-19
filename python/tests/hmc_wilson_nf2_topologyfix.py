from pybridge import *

setup(lattice_size=[4,4,4,8], grid_size=[], verbose_level="Detailed")

nvol=Common.Nvol
ndim=Common.Ndim

#-- gauge part

action_g=Action_G_Rectangle(beta=2.3,
                            c_plaq=3.648,
                            c_rect=-0.331)

#-- smearing part

smear_params=dict(rho_uniform=0.1, number_of_smearing=2)

proj=Projection_Stout_SU3()
smear=Smear_APE(proj, **smear_params)
dr_smear=Director_Smear(smear, **smear_params)

smearforce=ForceSmear_APE(proj, **smear_params)


#-- Nf=2 part

fopr_params=dict(gamma_matrix_type="Chiral",
                 hopping_parameter=0.5/(4.0-1.6),
                 boundary_condition=[1,1,1,-1])

fopr=Fopr_Wilson(**fopr_params)
force_fopr=Force_F_Wilson_Nf2(**fopr_params)

fopr_smear=Fopr_Smeared(fopr, dr_smear)
force_fopr_smear=Force_F_Smeared(force_fopr, smearforce, dr_smear)


solver_md=Solver_BiCGStab_L_Cmplx(fopr_smear,
                                  maximum_number_of_iteration=100,
                                  maximum_number_of_restart=40,
                                  convergence_criterion_squared=1.0e-28,
                                  use_initial_guess="false",
                                  Omega_tolerance=0.60,
                                  number_of_orthonormal_vectors=2,
                                  tolerance_for_DynamicSelection_of_L=0.1)
fprop_md=Fprop_Standard_lex(solver_md)

solver_h=Solver_BiCGStab_L_Cmplx(fopr_smear,
                                 maximum_number_of_iteration=100,
                                 maximum_number_of_restart=40,
                                 convergence_criterion_squared=1.0e-28,
                                 use_initial_guess="false",
                                 Omega_tolerance=0.60,
                                 number_of_orthonormal_vectors=2,
                                 tolerance_for_DynamicSelection_of_L=0.1)
fprop_h=Fprop_Standard_lex(solver_h)


#-- preconditioner part

fopr_prec_params=dict(gamma_matrix_type="Chiral",
                      hopping_parameter=0.5/(4.0-1.6),
                      twisted_mass=0.1,
                      boundary_condition=[1,1,1,-1])

fopr_prec=Fopr_Wilson_TwistedMass(**fopr_prec_params)
force_fopr_prec=Force_F_Wilson_TwistedMass_Nf2(**fopr_prec_params)

fopr_smear_prec=Fopr_Smeared(fopr_prec, dr_smear)
force_fopr_smear_prec=Force_F_Smeared(force_fopr_prec, smearforce, dr_smear)


solver_h_prec=Solver_BiCGStab_L_Cmplx(fopr_smear_prec,
                                      maximum_number_of_iteration=100,
                                      maximum_number_of_restart=40,
                                      convergence_criterion_squared=1.0e-28,
                                      use_initial_guess="false",
                                      Omega_tolerance=0.60,
                                      number_of_orthonormal_vectors=2,
                                      tolerance_for_DynamicSelection_of_L=0.1)
fprop_h_prec=Fprop_Standard_lex(solver_h_prec)


#-- fermion action

action_f=Action_F_Ratio_lex(
    fopr_smear_prec, force_fopr_smear_prec,
    fopr_smear, force_fopr_smear,
    fprop_h_prec,
    fprop_md, fprop_h)


#-- integrator

actions=[[action_f, action_g]]

directors=[dr_smear]

builder=Builder_Integrator(actions, directors,
                           integrator="Leapfrog",
                           #number_of_levels=1,
                           number_of_steps=[5],
                           order_of_exp_iP=8,
                           lambda_Omelyan=0.1931833275)
integ=builder.build()


#-- hmc
traj_number=200000

rand=RandomNumbers_Mseries(traj_number)

hmc=HMC_General(actions, directors, integ, rand,
                trajectory_length=0.1,
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
