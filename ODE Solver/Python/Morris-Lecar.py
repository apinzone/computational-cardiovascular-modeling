import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#Implement Morris Lecar Model

#Define Equations
def derivatives(t, y, params):
    #get state variables from initial array
    V, w = y
    #unpack parameters
    C, gCa, gK, gL, VCa, VK, VL, V1, V2, V3, V4, I_ext = params
    #calculate relevant components for derivative computation
    m_inf = 0.5 * (1 + np.tanh((V-V1)/V2))
    w_inf = 0.5 * (1 + np.tanh((V-V3)/V4))
    tau_w = 1 / np.cosh((V-V3)/(2*V4))
    #compute derivatives
    dVdt = (1/C) * (I_ext - gCa * m_inf * (V - VCa) - gK * w * (V - VK) - gL * (V - VL))
    dwdt = (w_inf - w) / tau_w
    #return solutions
    return [dVdt, dwdt]

#Set parameters
C = 20
gCa = 4.4 
gK = 8.0
gL = 2
VCa = 120
VK = -84
VL = -60
V1 = -1.2
V2 = 18
V3 = 2
V4 = 30
I_ext = 100
#Put parameters in one array while legibly defining beforehand
base_params = [C, gCa, gK, gL, VCa, VK, VL, V1, V2, V3, V4, I_ext]

#Set intiial conditions and time_span (2 separate arrays)
y0 = [0, 0] #V, w
t_span = [0, 500]

#Call solver and any additional arguments
solved = solve_ivp(derivatives, t_span, y0, args = (base_params,), method = 'BDF', dense_output = True)

#Plot 
t_plot = np.linspace(0, 500, 5000)
y_plot = solved.sol(t_plot)
V_plot = y_plot[0]
w_plot = y_plot[1]

plt.figure()
plt.plot(t_plot, V_plot)
plt.xlabel('Time (s)')
plt.ylabel('Voltage (mV)')

plt.figure()
plt.plot(V_plot, w_plot)
plt.show()
