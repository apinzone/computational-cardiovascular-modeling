import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

#derivatives function
def derivatives(t, y, p):
    #unpack params
    mu = p[0]
    beta = p[1]
    gamma = p[2]

    #unpack state variables
    S = y[0]
    I = y[1]
    R = y[2]

    #relevant calcs
    N = S + I + R

    #compute derivatives
    dSdt = mu * N - beta * (S * I)/ N - (mu * S)
    dIdt = beta * (S * I) / N - (gamma * I) - (mu * I)
    dRdt = (gamma * I) - (mu * R)

    #return solutions
    return [dSdt, dIdt, dRdt]

#parameters [mu, beta, gamma]
pars = [0.02, 0.5, 0.1]

#initial conditions
#y0 = [S, I, R]
y0 = [990, 10, 0]
t_span = [0, 200]

#Call solver
output = solve_ivp(derivatives, t_span, y0, args = (pars,), method = 'BDF', dense_output = True)
t_plot = np.linspace(0, 200, 2000)
y_plot = output.sol(t_plot)
S_plot = y_plot[0]
I_plot = y_plot[1]
R_plot = y_plot[2]

plt.plot(t_plot, S_plot, label = 'S')
plt.plot(t_plot, I_plot, label = 'I')
plt.plot(t_plot, R_plot, label = 'R')
plt.legend()
plt.xlabel('Time')
plt.show()