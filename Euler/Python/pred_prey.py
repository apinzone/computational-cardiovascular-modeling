import numpy as np
import matplotlib.pyplot as plt
import scipy 

#Define two equations
def dpdt (P, alpha, beta, Q):
    return alpha * P - beta * P * Q

def dqdt (P, delta, gamma, Q):
    return delta * P * Q - gamma * Q

#Set Parameters
alpha = 1 #prey birth rate
beta = 0.1 #predation rate
delta = 0.075 #predator growth rate from eating prey
gamma = 1.5 #predator death rate

#Set initial conditions 
P = 10 #starting prey population
Q = 5.0 #starting predator population
t = 0
t_end = 50
dt = 0.01

#create arrays
prey, pred, times = [P], [Q], [t]

#Solver loop
while t < t_end:
    dp = dpdt (P, alpha, beta, Q) * dt #don't forget to multiply time step
    dq = dqdt(P, delta, gamma, Q) * dt 
    P = P + dp #update prey state
    Q = Q + dq #update pred state
    t = t + dt #update time

    #update arrays
    prey.append(P)
    pred.append(Q)
    times.append(t)

plt.plot(times, prey, label = 'prey')
plt.plot(times, pred, label = 'predators')
plt.legend()
plt.xlabel('Time')
plt.ylabel('Populations')
plt.plot(prey, pred)
plt.show()
