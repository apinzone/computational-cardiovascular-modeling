import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#Define Derivative Function
def derivatives(t, y):
    #Get state variables from initial array
    V = y[0] #get V from y0 (passed into function)
    m = y[1] #get m from y0 (passed into function)
    h = y[2] #get h from y0
    n = y[3] #get n from y0

    #calculate currents and gates 
    INa = gNa * m**3 * h * (V - ENa)
    IK = gK * n**4 * (V - EK)
    IL = gL * (V- EL)
    alpha_m_v = 0.1 * (V + 40) / (1 - np.exp(-(V + 40) / 10))
    beta_m_v = 4.0 * np.exp(-(V +65)/18)
    alpha_h_v = 0.07 * np.exp(-(V +65)/20)
    beta_h_v = 1.0 / (1 + np.exp(-(V+35)/10))
    alpha_n_v = 0.01 * (V + 55) / (1 - np.exp(-(V + 55)/10))
    beta_n_v = 0.125 * np.exp(-(V + 65)/80)

    #Compute Derivatives
    dVdt = (1/Cm) * (I_ext - INa - IK - IL)
    dmdt = alpha_m_v * (1 - m) - beta_m_v * m
    dhdt = alpha_h_v * (1 - h) - beta_h_v  * h
    dndt = alpha_n_v * (1 - n) - beta_n_v * n
    #Update array all at once
    return [dVdt, dmdt, dhdt, dndt]

#Parameters (constants)
Cm = 1.0 #membrane capacitance
gNa = 120 #max sodium conductance
gK = 36 #max potassium conductance
gL = 0.3 #leak conductance
ENa = 50 #sodium reversal
EK = -77.0 #potassium reversal
EL = -54.4 #leak reversal
I_ext = 10 #Ext stimulus current

#Time Span 
t_span = [0, 50]

#Initial Conditions Array (y0)
y0 = [-65.0,0.05, 0.6, 0.32] # V, m ,h, n

#CALL SOLVER
solved = solve_ivp(derivatives, t_span, y0, method = 'BDF', dense_output = True)

#PLOT IT
t_plot = np.linspace(0, 50, 5000)
y_plot = solved.sol(t_plot)
plt.plot(t_plot, y_plot[0])
plt.xlabel('Time (s)')
plt.ylabel('Voltage (mV)')
plt.show()
