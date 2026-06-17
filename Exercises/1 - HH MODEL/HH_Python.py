import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#Rebuild Hodgkin-Huxley Model
#Initial Conditions and parameters
t_span = [0, 25] 
y0 = [-60, 0.3177, 0.0529, 0.5961]  #V, n, m, h
I_stim_range = np.linspace(0.05, 1.8, 15) #Range of 15 values from 0.05 to 1.8 to vary stimulus current
pulse_start = 2 #start pulse at 2ms 
pulse_end_range = [3, 3.25, 3.5, 3.75, 4, 5, 6, 7, 9, 12] #1ms, 1.25ms, 1.5ms, 1.75ms, 2ms, 3ms, 4ms, 5ms, 7 ms, 10 ms

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
    t_pulse_start = p[8] 
    t_pulse_end = p[9] 
    if t_pulse_start == 0 and t_pulse_end == 0: #Option for constant stimulus application 
        I_stim = p[7] 
    else:
        if t >= t_pulse_start and t <= t_pulse_end: #Option for pulse
            I_stim = p[7] 
        else:   
            I_stim = 0 

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

plt.figure() #Start figure 1 outside loop
#Solver Call across I_stim values 
for idx, I_stim in enumerate(I_stim_range):
    pars = [0.01, 55.17, -72.14, -49.42, 1.2, 0.36, 0.003, I_stim, 0, 0] 
           #C_m, E_Na, E_K, E_I, g_Na, g_K, g_leak, I_stim, blank pulse start and end times for constant stimulus
    output = solve_ivp(derivatives, t_span, y0, args = (pars,), method = 'BDF', dense_output = True) 
    t_plot = np.linspace(0, t_span[1], 200) 
    Y_values = output.sol(t_plot)
    V_array = Y_values[0]
    plt.subplot(3, 5, idx + 1)
    plt.plot(t_plot, V_array) #Figure 1. Plot V over time at diff stif current amplitudes 
    plt.xlabel('Time (ms)')
    plt.ylabel('Voltage (mv)')
    plt.ylim(-80, 60)
    plt.xlim(0, 25)
    plt.yticks(np.arange(-80, 61, 20))
    plt.xticks(np.arange(0, 26, 5))
    plt.title(f'I_stim = {I_stim:.3f}')
plt.tight_layout()

plt.figure() #Start figure 2 outside loop
#Solver Call across pulse durations
for idx, pulse_end_time in enumerate(pulse_end_range):
    pars_pulse = [0.01, 55.17, -72.14, -49.42, 1.2, 0.36, 0.003, 0.175, pulse_start, pulse_end_time] 
                 #C_m, E_Na, E_K, E_I, g_Na, g_K, g_leak, I_stim, pulse start time, pulse end time
    pulse_duration = pars_pulse[9] - pulse_start 
    output_pulse = solve_ivp(derivatives, t_span, y0, args = (pars_pulse,), method = 'BDF', dense_output = True) 
    t_plot_pulse = np.linspace(0, t_span[1], 200) 
    Y_values_pulse = output_pulse.sol(t_plot_pulse)
    V_array_pulse = Y_values_pulse[0]
    plt.subplot(2, 5, idx + 1)
    plt.plot(t_plot_pulse, V_array_pulse) #Figure 2, plot voltage over time at different pulse durations 
    plt.xlabel('Time (ms)')
    plt.ylabel('Voltage (mv)')
    plt.ylim(-80, 60)
    plt.xlim(0, 25)
    plt.yticks(np.arange(-80, 61, 20))
    plt.xticks(np.arange(0, 26, 5))
    plt.title(f'Pulse Duration = {pulse_duration:.1f}')
plt.tight_layout()

#Best parameters solver call for final figure
pars_best = [0.01, 55.17, -72.14, -49.42, 1.2, 0.36, 0.003, 0.175, 2, 3] #'Best' pulse duration of 1.5ms, I_stim of 0.175
    #1.C_m, 2.E_Na, 3.E_K, 4.E_I, 5.g_Na, 6.g_K, 7.g_leak, 8.stim current amp, 9. pulse start time, 10. pulse end time
output_best = solve_ivp(derivatives, t_span, y0, args = (pars_best,), method = 'BDF', dense_output = True) 
t_plot_best = np.linspace(0, t_span[1], 200) 
Y_values_best = output_best.sol(t_plot_best)
V_array_best = Y_values_best[0]

#Figure 3, plot voltage over time for 1 AP with 'Best' parameters
plt.figure()
plt.plot(t_plot_best, V_array_best) 
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (mv)')
plt.ylim(-80, 60)
plt.xlim(0, 25)
plt.yticks(np.arange(-80, 61, 20))
plt.xticks(np.arange(0, 26, 5))
plt.title('Action Potential: I stim = 0.175, Pulse Duration = 1 ms')

plt.show() #print all plots