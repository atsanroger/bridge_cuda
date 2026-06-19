from pybridge import *

setup(lattice_size=[4,4,4,8], grid_size=[], verbose_level="Detailed")

nvol=Common.Nvol
ndim=Common.Ndim

#-- gauge part

action_g=Action_G_Rectangle(beta=2.3,
                            c_plaq=3.648,
                            c_rect=-0.331)

#-- integrator

actions=[[action_g]]

directors=[]

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
