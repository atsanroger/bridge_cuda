from pybridge import *

setup(lattice_size=[4,4,4,8], grid_size=[], verbose_level="Detailed")

nvol=Common.Nvol
ndim=Common.Ndim

u=Field_G(nvol, ndim)

gconf=GaugeConfig("Text")
gconf.read_file(u, "conf_04040408.txt")


fopr=Fopr_Wilson(
          gamma_matrix_type="Dirac",
          hopping_parameter=0.5/(4.0-1.6),
          boundary_condition=[1,1,1,1])
fopr.set_config(u)

source=Source("Local",
              source_position=[0,0,0,0],
)


nshift=4

fprop=Fprop_Wilson_Shift(fopr,
                         number_of_shifts=nshift,
                         shifted_mass_difference=[0.10, 0.20, 0.30, 0.40],
                         maximum_number_of_iteration=4000,
                         convergence_criterion_squared=1.0e-24)

nc=Common.Nc
nd=Common.Nd

#for ispin in range(nd):
#    for icolor in range(nc):
for ispin in range(1):
    for icolor in range(1):
        b=Field_F(nvol)
        source.set(b, icolor=icolor, ispin=ispin)

        xs=[]
        for i in range(nshift):
          x=b.clone()
          xs.append(x)

        result=fprop.invert_D(xs, b)

