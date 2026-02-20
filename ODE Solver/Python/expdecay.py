import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#define equations
def derivatives(t, y):
    #Get state variables
    v, y_var = y
    #Compute derivs
    dvdt = -g - k * v
    dydt = v 
    #update array with derivatives
    return [dvdt, dydt]

#Define parameters
g = 9.8 #gravity constant
k = 0.5 #drag coefficient

#Define initial conditions and timespan
time_span = [0,10]
y0 = [0 , 100] #initial values for v, y

#Call solver 
output = solve_ivp(derivatives, time_span, y0, method = 'BDF', dense_output = True)

#plot 
t_plot = np.linspace(0, 10, 100)
y_plot = output.sol(t_plot)
v_plot = y_plot[0]
y_var_plot = y_plot[1]
plt.figure()
plt.plot(t_plot, v_plot)
plt.figure()
plt.plot(t_plot, y_var_plot)
plt.show()