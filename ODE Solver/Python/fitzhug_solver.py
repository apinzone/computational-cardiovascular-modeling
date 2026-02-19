import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#Define Equations
def derivatives(y, t):
    #Extract state variables
    V = y[0]
    W = y[1]

    #compute derivatives
    dVdt = V - (V**3/3) - W + I
    dWdt = (V + 0.7 - 0.8 * W)/12.5

    #Pass solutions as single array
    return [dVdt, dWdt]

#Set parameters
I = 0.5 #Ext stimulus current
#Set initial conditions in a single array (y0)
y0 = [-1.0, 0.0] #V, W

#Timespan
t_span = [0, 100]

#Call Solver
solutions = solve_ivp(derivatives, y0, t_span, method = 'BDF', dense_output = True)

#Plot 
time_array = np.linspace(0, 100, 1000)
y_array = sol.solutions(t_plot)
voltage_plot = y_array[0]
plt.plot(time_array,voltage_plot)
plt.xlabel('time (s)')
plt.ylabel('Voltage (mV)')
plt.show()