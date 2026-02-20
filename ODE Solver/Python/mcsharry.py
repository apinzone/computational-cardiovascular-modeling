import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#McSharry Dynamical Model
#Define Equations 
def derivatives(t, y, omega):
    #obtain state variables
    x, y_var, z = y
    #define relevant parameters 
    alpha = 1 - np.sqrt(x**2 + y_var**2)
    # omega = (2 * np.pi) #constant at 60 bpm, no peak detection 
    #current angle on circle
    theta = np.arctan2(y_var, x)
    #parameters for P, Q, R, S, T
    a_sub_i = np.array([1.2, -5.0, 30.0, -7.5, 0.75]) #P, Q, R, S, T
    b_sub_i = np.array([0.25, 0.1, 0.1, 0.1, 0.4]) #P, Q, R, S, T
    theta_i = np.array([-np.pi/3, -np.pi/12, 0, np.pi/12, np.pi/2])
    z0 = 0.005 #baseline value

    #compute sum contributions of PQRST for Z dot 
    sum_term = 0
    for i in range(5):
        delta_theta = np.mod(theta - theta_i[i] + np.pi, 2*np.pi) - np.pi
        sum_term += a_sub_i[i] * delta_theta * np.exp(-(delta_theta**2/(2*b_sub_i[i]**2)))
                                                
    #compute derivatives
    x_dot = (alpha * x) - (omega * y_var)
    y_dot = (alpha * y_var) + (omega * x)
    z_dot =  -sum_term - (z - z0)

    #return derivatives as a single array
    return[x_dot, y_dot, z_dot]

#

#Define initial conditions array and t_span
y0 = [1, 0, 0] #x, y_var, z
t_span = [0, 15]

#Varied parameters omega
#create range of omega values 
omega_values = np.linspace(np.pi, np.pi * 4, 20)

#array for storage
z_array = [] 

#Loop over each omega value 
for omega in omega_values:
    #run simulation
    solutions = solve_ivp(derivatives, t_span, y0, 
    args = (omega,), method = 'BDF', dense_output = True)
    #z trace
    z = solutions.y[2, :]
    #measure QRS amplitude
    qrs_amp = np.max(z) - np.min(z)
    #store appended array
    z_array.append(qrs_amp)

# #Call Solver 
# solutions = solve_ivp(derivatives, t_span, y0, method = 'BDF', dense_output = True)
# t_plot = np.linspace(1.5, 300, 30000)
# y_plot = solutions.sol(t_plot)
# x_plot = y_plot[0]
# y_var_plot = y_plot[1]
# z_plot = y_plot[2]

plt.figure()
plt.plot(omega_values, z_array, 'o-')
plt.xlabel('Omega (rad/s)')
plt.ylabel('QRS Amplitude (mV)')
plt.title('Heart Rate vs QRS Amplitude')
plt.grid(True)
plt.show()

# #Plot Generated ECG (Z)
# plt.figure()
# plt.plot(t_plot, z_plot)
# plt.xlabel('Time (s)')
# plt.ylabel('ECG (mV)')
# plt.show()