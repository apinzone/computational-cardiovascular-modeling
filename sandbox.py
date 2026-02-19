import numpy as np
import matplotlib.pyplot as plt
import scipy 

# #Define two equations
# def dpdt (P, alpha, beta, Q):
#     return alpha * P - beta * P * Q

# def dqdt (P, delta, gamma, Q):
#     return delta * P * Q - gamma * Q

# #Set Parameters
# alpha = 1 #prey birth rate
# beta = 0.1 #predation rate
# delta = 0.075 #predator growth rate from eating prey
# gamma = 1.5 #predator death rate

# #Set initial conditions 
# P = 10 #starting prey population
# Q = 5.0 #starting predator population
# t = 0
# t_end = 50
# dt = 0.01

# #create arrays
# prey, pred, times = [P], [Q], [t]

# #Solver loop
# while t < t_end:
#     dp = dpdt (P, alpha, beta, Q) * dt #don't forget to multiply time step
#     dq = dqdt(P, delta, gamma, Q) * dt 
#     P = P + dp #update prey state
#     Q = Q + dq #update pred state
#     t = t + dt #update time

#     #update arrays
#     prey.append(P)
#     pred.append(Q)
#     times.append(t)

# # plt.plot(times, prey, label = 'prey')
# # plt.plot(times, pred, label = 'predators')
# # plt.legend()
# # plt.xlabel('Time')
# # plt.ylabel('Populations')
# plt.plot(prey, pred)
# # plt.show()

#Van der Pol Oscillator Example
#Define Equations
def dxdt(y):
    return y

def dydt(x, y, mu):
    return mu * (1-x**2) * y - x

#Set Parameters
mu = 1.0

#Set initial conditions
x = 0.1
y = 0
t = 0
dt = 0.01
t_end = 30

#Create arrays for updating states
time_array, x_array, y_array = [t], [x], [y]
#Solver Loop call
while t<t_end:
    #compute derivatives
    dx = dxdt(y) * dt
    dy = dydt(x, y, mu) * dt
    #update state variables
    x = x + dx
    y = y + dy
    t = t + dt
    #update arrays
    time_array.append(t)
    x_array.append(x)
    y_array.append(y)

#Plot Results

plt.plot(time_array, x_array)
plt.plot(time_array, y_array)
plt.plot(x_array, y_array)
plt.show()

