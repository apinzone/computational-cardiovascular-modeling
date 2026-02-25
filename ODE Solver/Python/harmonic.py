import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#Define Equations
def derivatives(t, y, p):
    #unpack params
    k, m, c = p
    #Initial conditions
    x, v = y
    #compute derivatives
    dxdt = v
    dvdt = -(k/m) * x - (c/m) * v
    #return derivs as a single array
    return [dxdt, dvdt]

#params [k, m, c]
pars = [10, 1, 0.5]

#Initial conditions and timespan
y0 = [1, 0] #[x_0, v_0]
t_span = [0, 20]

#Solver call 
output = solve_ivp(derivatives, t_span, y0, args = (pars,), method = 'BDF', dense_output = True)

#plot
t_plot = np.linspace(0, 20, 200)
y_plot = output.sol(t_plot)
plt.plot(t_plot, y_plot[0]) #x over time
plt.plot(t_plot, y_plot[1]) #v over time
plt.show()