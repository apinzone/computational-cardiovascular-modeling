import random
from deap import base, creator, tools, algorithms

#Define Fitness (Minimization Problem)
creator.create("FitnessMin", base.Fitness, weights = (-1.0,))
creator.create("Individual", list, fitness=creator.FitnessMin)

#Toolbox Setup
toolbox = base.Toolbox()

#Each individual = one x value between -10 and 10
toolbox.register("attr_float", random.uniform, -10,10)
toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_float, n=1)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

#Fitness Function J(theta)
def evaluate(individual):
        x = individual[0]
        return (x-3)**2, #Tuple

#Register GA Operators
toolbox.register("evaluate", evaluate)
toolbox.register("mate", tools.cxBlend, alpha = 0.5)
toolbox.register("mutate", tools.mutGaussian, mu = 0, sigma = 1, indpb = 1)
toolbox.register("select", tools.selTournament, tournsize = 3)

#Run GA
population = toolbox.population(n = 50) #50 individuals
result, log = algorithms.eaSimple(
        population, toolbox,
        cxpb = 0.5,
        mutpb = 0.2,
        ngen = 20,
        verbose = True
)

# Best solution found
best = tools.selBest(result, k=1)[0]
print(f"\nBest x found: {best[0]:.4f}")
print(f"f(x) = {evaluate(best)[0]:.6f}")
print(f"True minimum: x = 3, f(x) = 0")

