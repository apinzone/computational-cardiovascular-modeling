import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#Implement Pendulum Swing Model with Damping

#Create Derivatives Function
def derivatives(t, y, p):
    #unpack params 
    g = p[0]
    L = p[1]
    b = p[2]

    #unpack state variables
    theta = y[0]
    omega = y[1]

    #compute derivatives
    theta_dot = omega #theta 
    omega_dot = -(g/L) * np.sin(theta) - b * omega #omega
    return [theta_dot, omega_dot]

#Set parameters [g , L , b]
pars = [9.8, 1.0, 0.5]

#Set initial conditions and time span [theta, omega]
y0 = [np.pi/2 , 0]
t_span = [0,20]

#Call Solver 
output = solve_ivp(derivatives, t_span, y0, args = (pars,), method = 'BDF', dense_output = True)

#Plot results 
t_plot = np.linspace(0, 20, 200)
y_plot = output.sol(t_plot)
theta_plot = y_plot[0]
omega_plot = y_plot[1]

#plot angle and velocity over time
plt.plot(t_plot, theta_plot, label = 'theta')
plt.plot(t_plot, omega_plot, label = 'omega')
plt.legend()
plt.xlabel('Time (s)')
plt.show()
