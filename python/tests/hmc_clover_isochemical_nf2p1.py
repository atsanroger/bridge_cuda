from pybridge import *
import math  # for sqrt()

setup(lattice_size=[4,4,4,8], grid_size=[],
      random_number_type="Mseries",
      random_number_seed=1234567,
      #      verbose_level="Detailed")
      verbose_level="General")

nvol=Common.Nvol
ndim=Common.Ndim


solver_param_base=dict(#maximum_number_of_iteration=100,
                       maximum_number_of_restart=40,
                       #convergence_criterion_squared=1.0e-28,
                       use_initial_guess="false",
                       Omega_tolerance=0.60,
                       number_of_orthonormal_vectors=2,
                       tolerance_for_DynamicSelection_of_L=0.1)

#-- gauge part

action_g=Action_G_Rectangle(beta=2.3,
                            c_plaq=3.648,
                            c_rect=-0.331)

#-- smear part

smear_params=dict(rho_uniform=0.1, number_of_smearing=2)

proj=Projection_Stout_SU3()
smear=Smear_APE(proj, **smear_params)
dr_smear=Director_Smear(smear, **smear_params)

smearforce=ForceSmear_APE(proj, **smear_params)


#-- Nf=2 part

fopr_nf2_params=dict(type="Clover_Isochemical",
                     gamma_matrix_type="Dirac",
                     hopping_parameter=0.13,
                     clover_coefficient=1.0,
                     isospin_chemical_potential=0.3,
                     boundary_condition=[1,1,1,-1])

fopr_nf2=Fopr_Clover_Isochemical(**fopr_nf2_params)
force_fopr_nf2=Force_F_Clover_Nf2_Isochemical(**fopr_nf2_params)

fopr_smear_nf2=Fopr_Smeared(fopr_nf2, dr_smear)
force_fopr_smear_nf2=Force_F_Smeared(force_fopr_nf2, smearforce, dr_smear)

solver_nf2_md=Solver(type="BiCGStab_L_Cmplx",
                     op=fopr_smear_nf2,
                     maximum_number_of_iteration=100,
                     convergence_criterion_squared=1.0e-28,
                     **solver_param_base)
fprop_nf2_md=Fprop_Standard_lex(solver_nf2_md)

solver_nf2_h=Solver(type="BiCGStab_L_Cmplx",
                    op=fopr_smear_nf2,
                    maximum_number_of_iteration=100,
                    convergence_criterion_squared=1.0e-28,
                    **solver_param_base)
fprop_nf2_h=Fprop_Standard_lex(solver_nf2_h)


#-- Nf=2 part preconditioner

fopr_nf2_prec_params=dict(type="Clover_Isochemical",
                          gamma_matrix_type="Dirac",
                          hopping_parameter=0.11,
                          clover_coefficient=1.0,
                          isospin_chemical_potential=0.3,
                          boundary_condition=[1,1,1,-1])

fopr_nf2_prec=Fopr_Clover_Isochemical(**fopr_nf2_prec_params)
force_fopr_nf2_prec=Force_F_Clover_Nf2_Isochemical(**fopr_nf2_prec_params)

fopr_smear_nf2_prec=Fopr_Smeared(fopr_nf2_prec, dr_smear)
force_fopr_smear_nf2_prec=Force_F_Smeared(force_fopr_nf2_prec, smearforce, dr_smear)

solver_nf2_md_prec=Solver(type="BiCGStab_L_Cmplx",
                          op=fopr_smear_nf2_prec,
                          maximum_number_of_iteration=100,
                          convergence_criterion_squared=1.0e-28,
                          **solver_param_base)
fprop_nf2_md_prec=Fprop_Standard_lex(solver_nf2_md_prec)

solver_nf2_h_prec=Solver(type="BiCGStab_L_Cmplx",
                         op=fopr_smear_nf2_prec,
                         maximum_number_of_iteration=100,
                         convergence_criterion_squared=1.0e-28,
                         **solver_param_base)
fprop_nf2_h_prec=Fprop_Standard_lex(solver_nf2_h_prec)


#---- action

action_f_nf2_prec=Action_F_Standard_lex(
    fopr_smear_nf2_prec, force_fopr_smear_nf2_prec,
    fprop_nf2_md_prec, fprop_nf2_h_prec)

action_f_nf2_d=Action_F_Ratio_lex(
    fopr_smear_nf2_prec, force_fopr_smear_nf2_prec,
    fopr_smear_nf2, force_fopr_smear_nf2,
    fprop_nf2_h_prec,
    fprop_nf2_md, fprop_nf2_h)


#-- Nf=1 part

fopr_nf1_params=dict(type="Clover_Isochemical",
                     gamma_matrix_type="Dirac",
                     hopping_parameter=0.12,
                     clover_coefficient=1.0,
                     isospin_chemical_potential=0.3,
                     boundary_condition=[1,1,1,-1])

rational_params=dict(number_of_poles=16,
                     lower_bound=1.0e-2,
                     upper_bound=math.sqrt(10.0),
                     maximum_number_of_iteration=1000,
                     convergence_criterion_squared=1.0e-20)

fopr_nf1=Fopr_Clover_Isochemical(**fopr_nf1_params)
force_fopr_nf1=Force_F_Clover_Nf2_Isochemical(**fopr_nf1_params)

fopr_smear_nf1_langev=Fopr_Smeared(
    Fopr_Rational(fopr_nf1,
                  exponent_numerator=1,
                  exponent_denominator=4,
                  **rational_params),
    dr_smear)

fopr_smear_nf1_h=Fopr_Smeared(
    Fopr_Rational(fopr_nf1,
                  exponent_numerator=-1,
                  exponent_denominator=2,
                  **rational_params),
    dr_smear)
force_fopr_smear_nf1_md=Force_F_Smeared(
    Force_F_Rational(fopr_nf1, force_fopr_nf1,
                     exponent_numerator=-1,
                     exponent_denominator=2,
                     **rational_params),
    smearforce, dr_smear)
    
#---- action

action_f_nf1=Action_F_Rational(
    fopr_smear_nf1_langev,
    fopr_smear_nf1_h, force_fopr_smear_nf1_md)


#-- integrator

actions=[
    [action_f_nf2_d],
    [action_f_nf2_prec, action_f_nf1],
    [action_g]
]

directors=[dr_smear]

builder=Builder_Integrator(actions, directors,
                           integrator="Leapfrog",
                           #number_of_levels=3,
                           number_of_steps=[5,2,3],
                           order_of_exp_iP=8,
                           lambda_Omelyan=0.1931833275)
integ=builder.build()


#-- hmc
traj_number=200000

rand=RandomNumbers_Mseries(traj_number)

hmc=HMC_General(actions, directors, integ, rand,
                trajectory_length=0.1,
                Metropolis_test="true")

#-- measurement

ploop=PolyakovLoop()

rand_noise=RandomNumbers_Mseries(seed=200000)
noise=NoiseVector_Z2(rand_noise)

qsuscept=QuarkNumberSusceptibility_Wilson(fopr_smear_nf2, fprop_nf2_h, noise,
                                          number_of_noises=10)


#-- exec

u=Field_G(nvol,ndim)

gconf=GaugeConfig("Text")
gconf.read_file(u, "conf_04040408.txt")

gconfout=GaugeConfig("None") # for output


ntraj=1
ntraj_interval=100

print("HMC start: Ntraj = {}", ntraj)


result_ploop=ploop.measure(u)
print("Polyakov loop = {}".format(result_ploop))
    
fopr_smear_nf2.set_config(u)
result_qsuscept=qsuscept.measure()


for traj in range(ntraj):
    print("")
    print("traj = {}".format(traj))

    traj_id = traj_number + traj + 1

    r=hmc.update(u)

    if traj_id % ntraj_interval == 0:
        print("---- write config ----")
        gconfout.write_file(u, "conf_output.dat-{:06}".format(traj_id))

    result_ploop=ploop.measure(u)
    print("Polyakov loop = {}".format(result_ploop))
    
    fopr_smear_nf2.set_config(u)
    result_qsuscept=qsuscept.measure()
    

gconfout.write_file(u, "conf_output.dat")

print("---- done ----")
