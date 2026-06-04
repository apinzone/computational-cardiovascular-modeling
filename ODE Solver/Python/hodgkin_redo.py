import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#Derivatives Function
def derivatives(t, y, p):
    #Unpack parameters 
    I_app = p[0]
    gNa = p[1]
    gK = p[2] 
    gL = p[3] 
    ENa = p[4] 
    EK = p[5] 
    EL = p[6] 
    Cm = p[7] 
    #Unpack state variables
    V = y[0]
    m = y[1]
    h = y[2]
    n = y[3]
    #calculations for rate functions
    alpha_m_V = 0.1 * (V + 40) / (1 - np.exp(-(V + 40)/10)) 
    beta_m_V = 4.0 * np.exp(-(V+65)/18) 
    alpha_h_V = 0.07 * np.exp(-(V+65)/20)
    beta_h_V = 1 / (1 + np.exp(-(V+35)/10)) 
    alpha_n_V = 0.01 * (V + 55) / (1 - np.exp(-(V+55)/10)) 
    beta_n_V = 0.125 * np.exp(-(V + 65)/80) 
    #calculations for currents
    INa = gNa * m**3 * h * (V - ENa) 
    IK = gK * n**4 * (V - EK) 
    IL = gL * (V - EL) 
    #compute derivatives
    dVdt = (1/Cm) * (I_app - INa - IK - IL) 
    dmdt = alpha_m_V * (1 - m) - beta_m_V * m
    dhdt = alpha_h_V * (1 - h) - beta_h_V * h
    dndt = alpha_n_V * (1 - n) - beta_n_V * n
    derivs = [dVdt, dmdt, dhdt, dndt] 
    return derivs

#Initial Conditions
t_span = [0, 100] #Run for 100 seconds 
y0 = [-65, 0.05, 0.6, 0.32] #V_0, m_0, h_0, n_0 
I_app_range = np.linspace(0,20, 50)

#Blank Arrays for Vmax and min
V_max_array = []
V_min_array = [] 

#Parameters + Solver Call  
for i in I_app_range:
    #I_app, gNa, gK, gL, ENa, EK, EL, Cm
    pars = [i, 120, 36, 0.3, 50, -77, -54.4, 1.0] 
    output = solve_ivp(derivatives, t_span, y0, args = (pars,), method = 'BDF', dense_output = True)
    t_plot = np.linspace(50, 100, 100) #Simulate only last 50 s.
    Y_values = output.sol(t_plot)
    V_max_array.append(np.max(Y_values[0]))
    V_min_array.append(np.min(Y_values[0]))

#Plot bifurcation Diagram
plt.figure()
plt.plot(I_app_range, V_max_array, label = 'V Max', color = 'blue')
plt.plot(I_app_range, V_min_array, label = 'V Min', color = 'red')
plt.xlabel('I')
plt.ylabel('Voltage (mV)')
plt.legend()
plt.show()
