import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#define equations
def derivatives(t, y, k):
    #Get state variables
    v, y_var = y
    #Compute derivs
    dvdt = -g - k * v
    dydt = v 
    #update array with derivatives
    return [dvdt, dydt]

# #Define parameters
# g = 9.8 #gravity constant
# k = 0.5 #drag coefficient

# #Define initial conditions and timespan
# time_span = [0,10]
# y0 = [0 , 100] #initial values for v, y

# #Call solver 
# output = solve_ivp(derivatives, time_span, y0, method = 'BDF', dense_output = True)

#Implement parameter sweep
#Define fixed 
g = 9.8 #gravity constant

#Define initial conditions and timespan
time_span = [0,10]
y0 = [0 , 100] #initial values for v, y

#Create a varied parameter, in this case, k (the drag coefficient)
#1. define range of k value to test
k_values = np.linspace(0.1, 2.0, 20) #20 values ranging from 0.1 to 2.0

#2. storage in array
terminal_velocities = []

#3. loop over each value, calling solver each time 
for k in k_values:
    #run simulation with this k 
    solved = solve_ivp(derivatives, time_span, y0, args = (k,), method = 'BDF', dense_output = True)
    #final terminal velocity
    v_final = solved.y[0, -1] #last value of v 
    #store
    terminal_velocities.append(v_final)

#4. Plot sweep 
plt.plot(k_values, terminal_velocities, 'o-')
plt.xlabel('Drag coefficient k')
plt.ylabel('Terminal velocity (m/s)')
plt.title('Sensitivity of Terminal Velocity to Drag')
plt.grid(True)
plt.show()

# #plot 
# t_plot = np.linspace(0, 10, 100)
# y_plot = output.sol(t_plot)
# v_plot = y_plot[0]
# y_var_plot = y_plot[1]
# plt.figure()
# plt.plot(t_plot, v_plot)
# plt.figure()
# plt.plot(t_plot, y_var_plot)
# plt.show()