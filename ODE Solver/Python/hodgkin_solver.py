import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#Define Derivative Function
def derivatives(t, y, params):
    #Get state variables from initial array
    V, m, h, n = y  # get all state variables at once instead of one at a time

    #unpack parameters from array
    Cm, gNa, gK, gL, ENa, EK, EL, I_ext = params

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
base_params = [1.0 , 120, 36, 0.3, 50, -77, -54.4, 10]
                #Cm , gNa, gK, gL, ENa, EK, EL, I_ext

#Initial Conditions Array (y0) and Time Span
y0 = [-65.0,0.05, 0.6, 0.32] # V, m ,h, n
t_span = [0, 50]

#CALL SOLVER
solved = solve_ivp(derivatives, t_span, y0, args = (base_params,), 
method = 'BDF', dense_output = True)

#PLOT IT
t_plot = np.linspace(0, 50, 5000)
y_plot = solved.sol(t_plot)
plt.plot(t_plot, y_plot[0])
plt.xlabel('Time (s)')
plt.ylabel('Voltage (mV)')
plt.show()
