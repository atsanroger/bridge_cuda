import math
from pybridge import *

setup(lattice_size=[4,4,4,8], grid_size=[],
      random_number_type="Mseries",
      random_number_seed=1234567,
      verbose_level="Detailed")

xmin=1.0e-4
xmax=10.0

npoint=100

rational=Tools.Math_Rational(number_of_poles=16,
                             exponent_numerator=-1,
                             exponent_denominator=2,
                             lower_bound=1.0e-4,
                             upper_bound=10.0)

ex=-1.0/2.0

for i in range(npoint+1):
    x=xmin+(xmax-xmin)*i/npoint

    y1=rational.func(x)
    y2=math.pow(x, ex)

    print("x,rational,x^(n/d) = {:12.8f}  {:12.8f}  {:12.8f} "
          .format(x, y1, y2))

