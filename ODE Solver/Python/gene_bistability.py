import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#Derivative function
def derivatives(t, y, p):
    #Unpack parameters 
    alpha = p[0]
    beta = p[1]
    n = p[2]
    #State variables 
    A = y[0]
    B = y[1]
    #Compute derivatives
    dAdt = alpha / (1 + B**n) - A
    dBdt = alpha/ (1 + A**n) - B
    #Return Array of Derivatives
    return [dAdt, dBdt]

#Set up parameter sweeping and base params 
n_range = np.linspace(1, 5, 20) #at y0 [1.5, 1.5]
alpha_range = np.linspace(1, 4, 30) 
base_pars = [2.5, 1, 2] #alpha, beta, n no parameter sweep

#Initial conditions and timespan [A, B]
y0_B_high = [0.5, 3.0] 
y0_A_high = [3.0, 0.5] 
y0 = [1.5, 1.5] 
t_span = [0 , 50] 

#Call solver in specific cases 

#----- Normal Case at three different initial conditions ------
output_base = solve_ivp(derivatives, t_span, y0, args = (base_pars,), method = 'BDF', dense_output = True)
output_B_high = solve_ivp(derivatives, t_span, y0_B_high, args = (base_pars,), method = 'BDF', dense_output = True)
output_A_high = solve_ivp(derivatives, t_span, y0_A_high, args = (base_pars,), method = 'BDF', dense_output = True)

#Get plot variables
Shared_Time = np.linspace(0, 50, 500)
Y_plot_base = output_base.sol(Shared_Time)
Y_plot_B_high = output_B_high.sol(Shared_Time)
Y_plot_A_high = output_A_high.sol(Shared_Time)
A_base = Y_plot_base[0]
B_base = Y_plot_base[1]
A_B_high = Y_plot_B_high[0]
B_B_high = Y_plot_B_high[1]
A_A_high = Y_plot_A_high[0]
B_A_high = Y_plot_A_high[1]

#Create plots 
plt.figure()
plt.plot(Shared_Time, A_base, label = 'A at moderate', linestyle = 'dashed', color = 'blue')
plt.plot(Shared_Time, B_base, label = 'B at moderate', linestyle = 'dotted', color = 'red')
plt.plot(Shared_Time, A_B_high, label = 'A at high B', linestyle = 'dashed', color = 'green')
plt.plot(Shared_Time, B_B_high, label = 'B at high B', linestyle = 'dotted', color = 'purple')
plt.plot(Shared_Time, A_A_high, label = 'A at high A', linestyle = 'dashed', color = 'magenta')
plt.plot(Shared_Time, B_A_high, label = 'B at high A', linestyle = 'dotted', color = 'cyan')
plt.xlabel('Time')
plt.ylabel('Protein Concentration')
plt.legend()

#Hill Coefficient (n) Sweep
A_array = []
B_array = []
for i in n_range:
    pars = [2.5, 1, i]
    output = solve_ivp(derivatives, t_span, y0, args = (pars,), method = 'BDF', dense_output = True)
    Y_plot = output.sol(Shared_Time)
    Last_A = (Y_plot[0])[-1]
    Last_B = (Y_plot[1])[-1]
    A_array.append(Last_A)
    B_array.append(Last_B)

#B Coefficient Plots
plt.figure()
plt.plot(n_range, A_array, label = 'A', linestyle = 'solid', color = 'black')
plt.plot(n_range, B_array, label = 'B', linestyle = 'dashed', color = 'red')
plt.xlabel('n')
plt.ylabel('Protein Concentration')
plt.legend()


#Alpha Sweep
A_array_alpha_base = []
B_array_alpha_base = []
A_array_alpha_B = []
B_array_alpha_B = []
A_array_alpha_A = []
B_array_alpha_A = []
for i in alpha_range:
    pars = [i, 1, 2]
    #Base solver call 
    output_base = solve_ivp(derivatives, t_span, y0, args = (pars,), method = 'BDF', dense_output = True)
    #Base solver call B High
    output_B = solve_ivp(derivatives, t_span, y0_B_high, args = (pars,), method = 'BDF', dense_output = True)
    #Base solver call A High
    output_A = solve_ivp(derivatives, t_span, y0_A_high, args = (pars,), method = 'BDF', dense_output = True)
    #Get state variables and final values
    Y_plot_base = output_base.sol(Shared_Time)
    Y_plot_B = output_B.sol(Shared_Time)
    Y_plot_A = output_A.sol(Shared_Time)
    Last_A_base = (Y_plot_base[0])[-1]
    Last_B_base = (Y_plot_base[1])[-1]
    Last_A_B = (Y_plot_B[0])[-1]
    Last_B_B = (Y_plot_B[1])[-1]
    Last_A_A = (Y_plot_A[0])[-1]
    Last_B_A = (Y_plot_A[1])[-1]
    #Append appropriate arrays 
    A_array_alpha_base.append(Last_A_base)
    B_array_alpha_base.append(Last_B_base)
    A_array_alpha_B.append(Last_A_B)
    B_array_alpha_B.append(Last_B_B)
    A_array_alpha_A.append(Last_A_A)
    B_array_alpha_A.append(Last_B_A)

#Plot final values for A and B at alpha ranges
plt.figure()
plt.plot(alpha_range, A_array_alpha_base, label = 'A at moderate', linestyle = 'dashed', color = 'blue')
plt.plot(alpha_range, B_array_alpha_base, label = 'B at moderate', linestyle = 'dotted', color = 'red')
plt.plot(alpha_range, A_array_alpha_B, label = 'A at high B', linestyle = 'dashed', color = 'green')
plt.plot(alpha_range, B_array_alpha_B, label = 'B at high B', linestyle = 'dotted', color = 'purple')
plt.plot(alpha_range, A_array_alpha_A, label = 'A at high A', linestyle = 'dashed', color = 'magenta')
plt.plot(alpha_range, B_array_alpha_A, label = 'B at high A', linestyle = 'dotted', color = 'cyan')
plt.xlabel('Alpha')
plt.ylabel('Protein Concentration')
plt.legend()
plt.show()