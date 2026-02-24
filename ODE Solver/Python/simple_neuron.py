import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#derivatives function 
def derivatives(t, y, p, bm):
    #parameters based on condition
    I = 10 #constant
    if bm == 1:
        a = p[0]
        b = p[1]
        c = p[2]
        d = p[3]
    else:
        a = p[4]
        b = p[5]
        c = p[6]
        d = p[7]
    
    #state variables
    v = y[0]
    u = y[1]

    #compute derivatives
    dvdt = (0.04*v**2) + (5 * v) + 140 - u + I
    dudt = a * ((b*v) - u)

    #return derivatives
    return [dvdt, dudt]

#parameters [a, b, c, d (burst 1), a, b, c, d (burst 0)]
pars = [0.02, 0.2, -50, 2, 0.02, 0.2, -65, 8]

#initial conditions
y0 = [-65, -13] #V, U
t_span = [0, 300]

#call solver 
output = solve_ivp(derivatives, t_span, y0, args = (pars, 0,), method = 'BDF', dense_output = True)

#PLOT
t_plot = np.linspace(0, 300, 3000)
y_plot = output.sol(t_plot)
v_plot = y_plot[0]
# u_plot = y_plot[1]
plt.plot(t_plot, v_plot)
# plt.plot(t_plot, u_plot)
plt.show()
    