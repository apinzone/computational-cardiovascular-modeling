import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#Define Initial Conditions and I Range 
t_span = [0, 200] 
y0 = [-1, -0.5] #V0, w0
I_range = np.linspace(-0.5, 1.5, 7)

#Derivatives Function 
def derivatives(t, y, p):
    #Unpack state variables
    V =  y[0] 
    w = y[1] 
    #Unpack pars
    a = p[0] 
    b = p[1]
    epsilon = p[2]
    I = p[3]
    #Compute derivatives 
    dVdt = V - V**3/3 - w + I 
    dwdt = epsilon * (V + a - b * w) 
    derivs = [dVdt, dwdt]
    return derivs 

#Blank array for bifurcation plots 
bifurcation_array_V_max = np.zeros(len(I_range))
bifurcation_array_V_min = np.zeros(len(I_range))

#Call solver with pars 
for i in I_range:
    pars = [0.7, 0.8, 0.08, i]
    output = solve_ivp(derivatives, t_span, y0, args = (pars,), method = 'BDF', dense_output = True)
    t_plot = np.linspace(100, 200, 1000) #Simulate only last 100 seconds to trim
    Y_values = output.sol(t_plot)
    V_max = np.max(Y_values[0])
    V_min = np.min(Y_values[0])
    bifurcation_array_V_max.append(V_max)
    bifurcation_array_V_min.append(V_min)

#Bifucation plot
plt.figure()
plt.plot(I_range, bifurcation_array_V_max, label = 'V Max', color = 'blue')
plt.plot(I_range, bifurcation_array_V_min, label = 'V Min', color = 'red')
plt.xlabel('I Values')
plt.ylabel('Voltage (mV)')
plt.legend()
plt.show()