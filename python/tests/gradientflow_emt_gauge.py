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


emt=EnergyMomentumTensor(c_plaq=5.0/3.0,
                         c_rect=-1.0/12.0,
                         max_momentum=1)

action=Action_G_Rectangle(c_plaq=1.0, c_rect=0.0, beta=float(Common.Nc))

gflow=GradientFlow(action,
                   order_of_RungeKutta=4,
                   step_size=0.02,
                   order_of_approx_for_exp_iP=8,
                   adaptive="false",
                   tolerance=1.0e-4,
                   safety_factor=0.9)

nstep=10
max_flow_time=0.2

tflow=0.0
for i in range(nstep):
    tflow,result=gflow.evolve(tflow,u)

    emt.set_field_strength(u)

    emt.measure_EMT(tflow)
    emt.measure_EMT_at_t(tflow)
    emt.measure_EMT_at_t_FT(tflow)
    emt.measure_EMT_at_x(tflow)
    emt.measure_EMT_at_x_FT(tflow)
    emt.measure_EMT_at_y(tflow)
    emt.measure_EMT_at_y_FT(tflow)
    emt.measure_EMT_at_z(tflow)
    emt.measure_EMT_at_z_FT(tflow)

    if tflow > max_flow_time:
        break
