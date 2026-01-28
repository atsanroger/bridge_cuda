import argparse

parser=argparse.ArgumentParser(description="test GradientFlow")

parser.add_argument("-v", "--verbose", help="increase verbosity",
                    action="count", default=0)
parser.add_argument("-q", "--quiet", help="decrease verbosity",
                    action="count", default=0)

parser.add_argument("--rk", help="order of Runge-Kutta",
                    type=int, choices=[1,2,3,4], default=3)

args=parser.parse_args()


rk_order=args.rk

vl = max(min(args.verbose-args.quiet+1, 3), 0)
vlevel = ["Crucial", "General", "Detailed", "Paranoiac"][vl]


from pybridge import *

setup(lattice_size=[4,4,4,8], grid_size=[], verbose_level=vlevel)

nvol=Common.Nvol
ndim=Common.Ndim

u=Field_G(nvol,ndim)

gconf=GaugeConfig("Text")
gconf.read_file(u, "conf_04040408.txt")


energy_density=EnergyDensity(c_plaq=1.0, c_rect=0.0)

action=Action_G_Rectangle(c_plaq=1.0, c_rect=0.0, beta=float(Common.Nc))

gflow=GradientFlow(action,
                   order_of_RungeKutta=rk_order,
                   step_size=0.01,
                   order_of_approx_for_exp_iP=8,
                   adaptive="false",
                   tolerance=1.0e-4,
                   safety_factor=0.9)

nstep=4
max_flow_time=1.0

tflow=0.0
for i in range(nstep):
    tflow,result=gflow.evolve(tflow,u)

    t2=tflow**2
    e_plaq=energy_density.E_plaq(u)
    e_clover=energy_density.E_clover(u)
    
    print("  (t, t^2 E_plaq, t^2 E_clover) = {:.8f} {:.16f} {:.16f}"
          .format(tflow, t2*e_plaq, t2*e_clover))

    if tflow > max_flow_time:
        break
