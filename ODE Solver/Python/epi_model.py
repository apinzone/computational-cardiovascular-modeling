import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#Epidemic Model
#Define derivative functions ALL IN ONE PLACE
def derivatives(t, y):
    #extract state variables
    S = y[0]
    I = y[1]
    R = y[2]
    #compute derivatives
    dSdt = -beta * S * I
    dIdt = beta * S * I - gamma * I
    dRdt = gamma * I
    #return all derivatives in a single array 
    return[dSdt,dIdt,dRdt]

#Parameters
beta = 0.5
gamma = 0.1 

#Initial conditions (one array for both variables) and timestep
y0 = [0.99, 0.01, 0.0]
t_span = [0, 100]

#Solver call
solved = solve_ivp(derivatives, t_span, y0, method = 'RK45', dense_output = True)

#plotting 
t_plot = np.linspace(0, 100, 1000)
y_plot = sol.solved(t_plot)
S_plot = y_plot[0]
I_plot = y_plot[1]
R_plot = y_plot[2]
plt.plot(t_plot, S_plot, label='Susceptible')
plt.plot(t_plot, I_plot, label='Infected')
plt.plot(t_plot, R_plot, label='Recovered')
plt.xlabel('Time (days)')
plt.ylabel('Population fraction')
plt.legend()
plt.show()