
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#Define Parameters
#create sweep for V_in
V_in_array = np.linspace(0.1, 1.6, 16)
#Initial conditions 
y0 = [10, 5] 
t_span = [0, 500]

#Derivatives function
def derivatives(t, y, p):
    #unpack pars and vary 
    V_in, k_1, k_p, k_m = p
    #state variables
    G, ATP = y
    #compute derivatives
    G_dot = V_in - k_1 * G * ATP
    ATP_dot = 2*k_1 * G * ATP - ((k_p * ATP)/(ATP+k_m))
    #return array
    return [G_dot, ATP_dot]

#Create arrays for bifurcation plot 
max_glucose = np.zeros(len(V_in_array))
min_glucose = np.zeros(len(V_in_array))
max_ATP = np.zeros(len(V_in_array))
min_ATP = np.zeros(len(V_in_array))

#solver loop
for idx, v in enumerate(V_in_array):
    #define parameters with V_in varied 
    #[V_in, K_1, k_p, K_m]
    pars = [v, 0.02, 6, 13]
    #call solver 
    solutions = solve_ivp(derivatives, t_span, y0, args = (pars,), method = 'BDF', dense_output = True)
    time_plot = np.linspace(0, 500, 5000)
    y_plot = solutions.sol(time_plot)
    time_20_idx = int(len(time_plot) * 0.8)
    max_glucose[idx] = np.max(y_plot[0, time_20_idx:])  
    min_glucose[idx] = np.min(y_plot[0, time_20_idx:])
    max_ATP[idx] = np.max(y_plot[1, time_20_idx:])
    min_ATP[idx] = np.min(y_plot[1, time_20_idx:]) 

plt.figure(figsize=(10, 8))

plt.subplot(2, 1, 1)
plt.plot(V_in_array, max_glucose, 'b.-', markersize=8, label='Max')
plt.plot(V_in_array, min_glucose, 'r.-', markersize=8, label='Min')
plt.xlabel('$V_{in}$', fontsize=12)
plt.ylabel('[Glucose]', fontsize=12)
plt.title('Bifurcation Diagram: Glucose vs $V_{in}$', fontsize=14)
plt.legend()
plt.grid(True, alpha=0.3)

plt.subplot(2, 1, 2)
plt.plot(V_in_array, max_ATP, 'b.-', markersize=8, label='Max')
plt.plot(V_in_array, min_ATP, 'r.-', markersize=8, label='Min')
plt.xlabel('$V_{in}$', fontsize=12)
plt.ylabel('[ATP]', fontsize=12)
plt.title('Bifurcation Diagram: ATP vs $V_{in}$', fontsize=14)
plt.legend()
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()
    # G_plot = y_plot[0]
    # ATP_plot = y_plot[1]

    # plt.figure()
    # plt.plot(time_plot, G_plot)
    # plt.plot(time_plot, ATP_plot)
    # plt.show()


