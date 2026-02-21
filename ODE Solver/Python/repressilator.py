import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#Derivatives Function
def derivatives(t, y, p):
    #unpack parameters
    alpha = p[0]
    alpha_0 = p[1]
    beta = p[2]
    n = p[3]

    #get state variables
    m1 = y[0]
    m2 = y[1]
    m3 = y[2]
    p1 = y[3]
    p2 = y[4]
    p3 = y[5]

    #compute derivatives
    dm1dt = -m1 + alpha/(1 + p3**n) + alpha_0
    dm2dt = -m2 + alpha/(1+p1**n) + alpha_0
    dm3dt = -m3 + alpha/(1+p2**n) + alpha_0
    dp1dt = -beta * (p1 - m1)
    dp2dt = -beta * (p2 - m2)
    dp3dt = -beta * (p3 - m3)

    #Return solved array
    return [dm1dt, dm2dt, dm3dt, dp1dt, dp2dt, dp3dt]

#Parameters 
alpha = 216
alpha_0 = 0
beta = 0.2
n = 2.0
base_params = [alpha, alpha_0, beta, n]

#Initial Conditions
t_span = [0, 100]
y0 = [0 , 0, 0, 2, 0, 0]
   #[m1, m2, m3, p1, p2, p3]

#call solver
output = solve_ivp(derivatives, t_span, y0, args = (base_params,), method = 'BDF', dense_output = True)

#Create plot variables
time_plot = np.linspace(0, 100, 1000)
y_plot = output.sol(time_plot)
m1_array = y_plot[0]
m2_array = y_plot[1]
m3_array = y_plot[2]

#Create plots
plt.plot(time_plot, m1_array, label = 'm1')
plt.plot(time_plot, m2_array, label = 'm2')
plt.plot(time_plot, m3_array, label = 'm3')
plt.xlabel('Time')
plt.ylabel('mRNA Concentration')
plt.legend()
plt.show()