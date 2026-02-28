import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#Define Equations
def derivatives(t, y, i):
    #Extract state variables
    V = y[0]
    W = y[1]

    #compute derivatives
    dVdt = V - (V**3/3) - W + i
    dWdt = (V + 0.7 - 0.8 * W)/12.5

    #Pass solutions as single array
    return [dVdt, dWdt]

#Ext stimulus current
#Set initial conditions in a single array (y0)
y0 = [-1.0, 0.0] #V, W
#Timespan
t_span = [0, 100]

#Vary parameters
I_ranges = np.linspace(0.3, 0.8, 5)
for i in I_ranges:
    #Call Solver
    solutions = solve_ivp(derivatives, t_span, y0, args = (i,), method ='BDF', dense_output = True)            
    #Plot 
    time_array = np.linspace(0, 100, 1000)
    y_array = solutions.sol(time_array)
    voltage_plot = y_array[0]
    w_plot = y_array[1]
    plt.figure()
    plt.plot(time_array,voltage_plot)
    plt.xlabel('time (s)')
    plt.ylabel('Voltage (mV)')
    plt.figure()
    plt.plot(voltage_plot, w_plot)
plt.show()


