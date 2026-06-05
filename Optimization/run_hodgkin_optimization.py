import random
from deap import base, creator, tools, algorithms
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#Derivatives Function
def derivatives(t, y, p):
    #Unpack parameters 
    I_app = p[0]
    gNa = p[1]
    gK = p[2] 
    gL = p[3] 
    ENa = p[4] 
    EK = p[5] 
    EL = p[6] 
    Cm = p[7] 
    #Unpack state variables
    V = y[0]
    m = y[1]
    h = y[2]
    n = y[3]
    #calculations for rate functions
    alpha_m_V = 0.1 * (V + 40) / (1 - np.exp(-(V + 40)/10)) 
    beta_m_V = 4.0 * np.exp(-(V+65)/18) 
    alpha_h_V = 0.07 * np.exp(-(V+65)/20)
    beta_h_V = 1 / (1 + np.exp(-(V+35)/10)) 
    alpha_n_V = 0.01 * (V + 55) / (1 - np.exp(-(V+55)/10)) 
    beta_n_V = 0.125 * np.exp(-(V + 65)/80) 
    #calculations for currents
    INa = gNa * m**3 * h * (V - ENa) 
    IK = gK * n**4 * (V - EK) 
    IL = gL * (V - EL) 
    #compute derivatives
    dVdt = (1/Cm) * (I_app - INa - IK - IL) 
    dmdt = alpha_m_V * (1 - m) - beta_m_V * m
    dhdt = alpha_h_V * (1 - h) - beta_h_V * h
    dndt = alpha_n_V * (1 - n) - beta_n_V * n
    derivs = [dVdt, dmdt, dhdt, dndt] 
    return derivs

#Initial Conditions for HH Model
t_span = [0, 100] #Run for 100 seconds 
y0 = [-65, 0.05, 0.6, 0.32] #V_0, m_0, h_0, n_0 

#Load experimental voltage trace 
data = np.loadtxt('voltage_trace.csv', delimiter=',')
data2 = np.loadtxt('voltage_trace2.csv', delimiter=',')
t_experimental = data[:, 0]
V_experimental = data[:, 1]
t_experimental2 = data2[:, 0]
V_experimental2 = data2[:, 1]

#Define Fitness (Minimization Problem)
creator.create("FitnessMin", base.Fitness, weights = (-1.0,))
creator.create("Individual", list, fitness=creator.FitnessMin)

#Toolbox Setup
toolbox = base.Toolbox()

#Physiologically Meaningful Params
def init_gNa(): return random.uniform(60,180)
def init_gK(): return random.uniform(18,54)
def init_gL(): return random.uniform(0.15, 0.45) 
def init_individual():
        return creator.Individual([init_gNa(), init_gK(), init_gL()])
toolbox.register("individual", init_individual)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

#Fitness Function J(theta)
def evaluate(individual):
        gNa, gK, gL = individual
        #Enforce bounds 
        if not (60 <= gNa <= 180): return 999999,
        if not (18 <= gK  <= 54):  return 999999,
        if not (0.15 <= gL <= 0.45): return 999999,
        #Protocol 1 - Standard I_app
        pars1 = [10, gNa, gK, gL, 50, -77, -54.4, 1.0]
        output1 = solve_ivp(derivatives, t_span, y0, args=(pars1,), method = 'BDF', dense_output = True)
        V_sim1 = output1.sol(t_experimental)[0]
        MSE1 = np.mean((V_sim1 - V_experimental)**2)

        #Protocol 2 - Higher I_app
        pars2 = [15, gNa, gK, gL, 50, -77, -54.4, 1.0]
        output2 = solve_ivp(derivatives, t_span, y0, args=(pars2,), method = 'BDF', dense_output = True)
        V_sim2 = output2.sol(t_experimental2)[0]
        MSE2 = np.mean((V_sim2 - V_experimental2)**2)

        MSE = MSE1 + MSE2 
        return MSE ,


#Register GA Operators
toolbox.register("evaluate", evaluate)
toolbox.register("mate", tools.cxBlend, alpha = 0.5)
toolbox.register("mutate", tools.mutGaussian, mu = 0, sigma = 0.5, indpb = 0.2)
toolbox.register("select", tools.selTournament, tournsize = 3)

#Run GA
population = toolbox.population(n = 200) #100 individuals
# Statistics tracking
stats = tools.Statistics(lambda ind: ind.fitness.values)
stats.register("min", np.min)
stats.register("avg", np.mean)

result, log = algorithms.eaSimple(
        population, toolbox,
        cxpb = 0.5,
        mutpb = 0.2,
        ngen = 200,
        stats = stats,
        verbose = True
)

# Extract and display best parameters
best = tools.selBest(result, k=1)[0]
gNa_found, gK_found, gL_found = best

print(f"\nTrue:  gNa=120,        gK=36,       gL=0.3")
print(f"Found: gNa={gNa_found:.2f}, gK={gK_found:.2f}, gL={gL_found:.4f}")
final_mse = evaluate(best)[0]
print(f"Final MSE: {final_mse:.4f}")

# Plot experimental vs recovered trace
pars_best = [10, gNa_found, gK_found, gL_found, 50, -77, -54.4, 1.0]
output_best = solve_ivp(derivatives, t_span, y0, 
                        args=(pars_best,), method='BDF', dense_output=True)
V_recovered = output_best.sol(t_experimental)[0]

plt.figure()
plt.plot(t_experimental, V_experimental, 'b-', label='Experimental', linewidth=2)
plt.plot(t_experimental, V_recovered, 'r--', label='GA Recovered', linewidth=2)
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (mV)')
plt.title('HH Parameter Recovery with DEAP')
plt.legend()
plt.show()
