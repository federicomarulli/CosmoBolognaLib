import numpy as np
import matplotlib.pyplot as plt
import CosmoBolognaLib as cbl
from CosmoBolognaLib import DoubleVector as dv

# Define the function to put on a grid
npoints = 200
xx = np.linspace(0., 4*np.pi, 200)
func = np.cos(xx)*np.exp(-0.1*xx)

# Define the b-spline
# number of breakpoints
nknots = 10

# b-spline order
order = 4

# Construct the b-spline object
interp_func = cbl.FuncGrid_Bspline(dv(xx), dv(func), nknots, order)

# Compute the interpolated values
xx2 = np.linspace(0., 4*np.pi, 50) 
func2 = interp_func.eval_func(dv(xx2))

# Plot results
plt.plot(xx, func, label="True function")
plt.plot(xx2, func2, label="b-spline")

plt.xlabel("x")
plt.ylabel("f(x)")

plt.legend(loc="best")
plt.show(block=False)
