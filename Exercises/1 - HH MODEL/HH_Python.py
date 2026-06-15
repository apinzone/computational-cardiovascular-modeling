import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#Rebuild Hodgkin-Huxley Model
#Initial Conditions and parameters
t_span = [0, 25] 
y0 = [-60, 0.3177, 0.0529, 0.5961]  #V, n, m, h
I_stim_range = np.linspace(0.05, 1.8, 15) #Range of 15 values from 0.05 to 1.8 to vary stimulus current

#Derivatives Function
def derivatives(t,y,p):
    #Unpack state variables 
    V = y[0]
    n = y[1]
    m = y[2]
    h = y[3] 
    #Unpack pars 
    C_m = p[0]
    E_Na = p[1] 
    E_K = p[2]
    E_leak = p[3]
    g_Na = p[4]
    g_K = p[5] 
    g_leak = p[6] 
    I_stim = p[7] 
    #Calculate rate constants 
    alpha_n = 0.01 * (V + 50) / (1 - np.exp(-(V + 50)/10)) 
    beta_n = 0.125 * np.exp((-(V+60))/80) 
    alpha_m = 0.1 * (V + 35) / (1 - np.exp(-(V + 35)/10)) 
    beta_m = 4 * np.exp(-0.0556*(V + 60)) 
    alpha_h = 0.07 * np.exp(-0.05 * (V + 60)) 
    beta_h = 1 / (1 + np.exp(-0.1 * (V + 30))) 
    #Calculate Currents
    I_Na = g_Na * m**3 * h * (V - E_Na) 
    I_K = g_K * n**4 * (V - E_K) 
    I_leak = g_leak * (V - E_leak) 
    #Compute Derivatives
    dVdt = (1/C_m) * (I_stim - I_Na - I_K - I_leak) 
    dndt = alpha_n * (1 - n) - beta_n * n  
    dmdt = alpha_m * (1 - m) - beta_m * m 
    dhdt = alpha_h * (1 - h) - beta_h * h 
    #Return derivatives
    return [dVdt, dndt, dmdt, dhdt]

plt.figure() #Start figure outside loop
#Solver Call across I_stim values 
for idx, I_stim in enumerate(I_stim_range):
    pars = [0.01, 55.17, -72.14, -49.42, 1.2, 0.36, 0.003, I_stim] 
           #C_m, E_Na, E_K, E_I, g_Na, g_K, g_leak, I_stim
    output = solve_ivp(derivatives, t_span, y0, args = (pars,), method = 'BDF', dense_output = True) 
    t_plot = np.linspace(0, t_span[1], 200) 
    Y_values = output.sol(t_plot)
    V_array = Y_values[0]
    plt.subplot(3, 5, idx + 1)
    plt.plot(t_plot, V_array) #Plot V over time
    plt.xlabel('Time (ms)')
    plt.ylabel('Voltage (mv)')
    plt.ylim(-80, 60)
    plt.xlim(0, 25)
    plt.yticks(np.arange(-80, 61, 20))
    plt.xticks(np.arange(0, 26, 5))
    plt.title(f'I_stim = {I_stim:.3f}')

plt.tight_layout()
plt.show() #print plot 




