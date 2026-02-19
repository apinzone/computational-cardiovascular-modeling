import numpy as np
import matplotlib.pyplot as plt
import scipy 

#CALCIUM TRANSIENT MODEL
#calcium derivative function
def dCadt (Ca, release, k_pump):
    return release - (k_pump * Ca)
#release derivative function
def dreleasedt (stimulus, release, tau):
    return (stimulus - release)/tau

#Parameters
k_pump = 0.5
stimulus = 1.0 
tau = 5.0

#Initial Conditions
Ca = 0.0
release = 0.0
t = 0
dt = 0.01 #small time step
t_end = 50

#Create arrays
times, calcium_array, release_array = [t], [Ca], [release]

#Solver lopp
while t < t_end:
    if t < 2: 
        stimulus = 1
    else: 
        stimulus = 0

    #compute derivatives
    dCa = dCadt(Ca, release, k_pump) * dt
    drelease = dreleasedt(stimulus, release, tau) * dt

    #Update state variables
    Ca = Ca + dCa
    release = release + drelease
    t = t + dt
    #Update arrays with derivatives at current point in time
    times.append(t)
    calcium_array.append(Ca)
    release_array.append(release)

plt.plot(times, calcium_array)
plt.xlabel('Time')
plt.ylabel('Calcium')
plt.plot(times, release_array)
plt.show()