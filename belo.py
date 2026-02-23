
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

#Belousov-Zhabotinsky Reaction

#Create derivative function
def derivatives(t, y, p):
    #parameters
    epsilon_1 = p[0]
    epsilon_2 = p[1]
    q = p[2]
    f = p[3]
    k_1 = p[4]
    k_2 = p[5]
    k_3 = p[6]
    k_4 = p[7]
    k_5 = p[8]
    #state variables
    x = y[0]
    y_var =y[1]
    z = y[2]
    w = y[3]
    v = y[4]
    #compute derivatives
    dxdt = (1/epsilon_1) * (q * y_var - x *y_var + x * (1-x))
    dydt = (1/epsilon_2) * (-q * y_var - x *y_var + f * z)
    dzdt = x -z 
    dwdt = -k_1 * w + k_2 * x - k_3 * w * x
    dvdt = k_4 * w * y_var - k_5 * v
    #return solution
    return [dxdt, dydt, dzdt, dwdt, dvdt]

#parameters [epsilon_1, epsilon_2, q, f, k1, k2, k3, k4, k5]
pars = [0.04, 1.0, 0.002, 1.0, 0.1, 0.2, 0.05, 0.3, 0.15]

#initial conditions [x, y_var, z, w, v]
y0 = [0.5, 0.1, 0.2, 0, 0]
t_span = [0, 500]

#call solver
output = solve_ivp(derivatives, t_span, y0, args = (pars,), method = 'BDF', dense_output = True)
t_plot = np.linspace(0, 500, 5000)
y_plot = output.sol(t_plot)
x_plot = y_plot[0]
y_var_plot = y_plot[1]
z_plot = y_plot[2]
w_plot = y_plot[3]
v_plot = y_plot[4]

plt.plot(t_plot, x_plot)
plt.plot(t_plot, y_var_plot)
plt.plot(t_plot, z_plot)
plt.plot(t_plot, w_plot)
plt.plot(t_plot, v_plot)
plt.show()