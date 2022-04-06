from multiprocessing import Pool
import pygad
from multiprocessing import Pool
import numpy as np
import numpy
last_fitness = 0
def f(x):
    return x*x

def run_wrap(pool):
    print (pool.map(f, range(10))  )        # prints "[0, 1, 4,..., 81]"

def run():
    pool = Pool(processes=4)     
    run_wrap(pool)# start 4 worker processes
    return 0

class PooledGA(pygad.GA):
    def __init__(self,*args,**kwargs):
        self.pool = Pool(processes=2)     
        super().__init__(*args,**kwargs)


    def cal_pop_fitness(self):
        pop_fitness = self.pool.map(fitness_wrapper, self.population)
        pop_fitness = np.array(pop_fitness)
        return pop_fitness



def fitness_func(solution, solution_idx):
    # Calculating the fitness value of each solution in the current population.
    # The fitness function calulates the sum of products between each input and its corresponding weight.
    function_inputs = [4,-2,3.5,5,-11,-4.7]
    output = np.sum(solution*function_inputs)
    desired_output = 44
    fitness = 1.0 / np.abs(output - desired_output)
    return fitness


def fitness_wrapper(solution):
    return fitness_func(solution, 0)

def callback_generation(ga_instance):
    global last_fitness
#     print("Generation = {generation}".format(generation=ga_instance.generations_completed))
#     print("Fitness    = {fitness}".format(fitness=ga_instance.best_solution()[1]))
#     print("Change     = {change}".format(change=ga_instance.best_solution()[1] - last_fitness))
    last_fitness = ga_instance.best_solution()[1]

def run_pga():
    function_inputs = [4,-2,3.5,5,-11,-4.7] # Function inputs.
    desired_output = 44 # Function output.



    fitness_function = fitness_func

    num_generations = 100 # Number of generations.
    num_parents_mating = 7 # Number of solutions to be selected as parents in the mating pool.

    # To prepare the initial population, there are 2 ways:
    # 1) Prepare it yourself and pass it to the initial_population parameter. This way is useful when the user wants to start the genetic algorithm with a custom initial population.
    # 2) Assign valid integer values to the sol_per_pop and num_genes parameters. If the initial_population parameter exists, then the sol_per_pop and num_genes parameters are useless.
    sol_per_pop = 5000 # Number of solutions in the population.
    num_genes = len(function_inputs)

    last_fitness = 0

    pool_ga = PooledGA(num_generations=num_generations,
                           num_parents_mating=num_parents_mating, 
                           fitness_func=fitness_function,
                           sol_per_pop=sol_per_pop, 
                           num_genes=num_genes,
                           on_generation=callback_generation)