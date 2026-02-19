import numpy as np
import matplotlib.pyplot as plt
import scipy 

# def dVdt(V):
#     return -V
# dt = 0.01
# t = 0
# V = 1.0
# t_end = 5.0

# times, voltages = [t], [V]

# while t<t_end:
#     V = V + dVdt(V) * dt
#     t = t + dt
#     times.append(t)
#     voltages.append(V)

# plt.plot(times, voltages)
# plt.xlabel('Time')
# plt.ylabel('Voltage')


# def dVdt(V, W, I):
#     return V - (V**3/3) - W + I
# def dWdt(W, V):
#     return (V + 0.7 - (0.8 *W))/12.5
# dt = 0.01
# I = 0.5
# V = -1.0
# W = 0.0
# t = 0
# t_end = 100

# times, voltages, recovery = [t], [V], [W]

# while t<t_end:
#     dV = dVdt(V, W, I) * dt
#     dW = dWdt(W, V) * dt 
#     V = V + dV
#     W = W + dW
#     t = t + dt
#     times.append(t)
#     voltages.append(V)
#     recovery.append(W)


#CALCIUM TRANSIENT MODEL
# def dCadt (Ca, release, k_pump):
#     return release - (k_pump * Ca)

# def dreleasedt (stimulus, release, tau):
#     return (stimulus - release)/tau

# k_pump = 0.5
# stimulus = 1.0 
# tau = 5.0
# Ca = 0.0
# release = 0.0
# t = 0
# dt = 0.01
# t_end = 50

# times, calcium_array, release_array = [t], [Ca], [release]

# while t < t_end:
#     if t < 2: 
#         stimulus = 1
#     else: 
#         stimulus = 0
#     dCa = dCadt(Ca, release, k_pump) * dt
#     drelease = dreleasedt(stimulus, release, tau) * dt
#     Ca = Ca + dCa
#     release = release + drelease
#     t = t + dt
#     times.append(t)
#     calcium_array.append(Ca)
#     release_array.append(release)

# plt.plot(times, calcium_array)
# plt.xlabel('Time')
# plt.ylabel('Calcium')
# plt.plot(times, release_array)
# plt.show()

#Implement Hodgkin-Huxley
#Define Equations 
def dVdt (V, gNa, gK, gL, m, h, n, ENa, EK, EL, Cm, I_ext):
    INa = gNa * m**3 * h * (V - ENa)
    IK = gK * n**4 * (V - EK)
    IL = gL * (V- EL)
    return (1/Cm) * (I_ext - INa - IK - IL)

def dmdt (m, V):
    alpha_m_v = 0.1 * (V + 40) / (1 - np.exp(-(V + 40) / 10))
    beta_m_v = 4.0 * np.exp(-(V +65)/18)
    return alpha_m_v * (1 - m) - beta_m_v * m

def dhdt (h, V):
    alpha_h_v = 0.07 * np.exp(-(V +65)/20)
    beta_h_v = 1.0 / (1 + np.exp(-(V+35)/10))
    return alpha_h_v * (1 - h) - beta_h_v  * h

def dndt (n, V):
    alpha_n_v = 0.01 * (V + 55) / (1 - np.exp(-(V + 55)/10))
    beta_n_v = 0.125 * np.exp(-(V + 65)/80)
    return alpha_n_v * (1 - n) - beta_n_v * n

#Define Parameters
Cm = 1.0 #membrane capacitance
gNa = 120 #max sodium conductance
gK = 36 #max potassium conductance
gL = 0.3 #leak conductance
ENa = 50 #sodium reversal
EK = -77.0 #potassium reversal
EL = -54.4 #leak reversal
I_ext = 10 #Ext stimulus current

#Set initial Conditions
V = -65.0
m = 0.05
h = 0.6
n = 0.32
t = 0
t_end = 50
dt = 0.01

times, voltage, m_array, h_array, n_array = [t], [V], [m], [h], [n]

#Solver Loop
while t < t_end:
    #solve each equation
    dv = dVdt (V, gNa, gK, gL, m, h, n, ENa, EK, EL, Cm, I_ext) * dt
    dm = dmdt (m, V) * dt
    dh = dhdt (h, V) * dt
    dn = dndt (n, V) * dt

    #Apply state variables to integrate forward in time
    V = V + dv
    m = m + dm
    h = h + dh
    n = n + dn
    t = t + dt

    #Append arrays
    times.append(t)
    voltage.append(V)
    m_array.append(m)
    h_array.append(h)
    n_array.append(n)

#Plot Voltage
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

ax1.plot(times, voltage, label='V')
ax1.set_ylabel('Voltage (mV)')
ax1.legend()

ax2.plot(times, m_array, label='m')
ax2.plot(times, h_array, label='h')
ax2.plot(times, n_array, label='n')
ax2.set_ylabel('Gate probability')
ax2.set_xlabel('Time')
ax2.legend()
plt.show()