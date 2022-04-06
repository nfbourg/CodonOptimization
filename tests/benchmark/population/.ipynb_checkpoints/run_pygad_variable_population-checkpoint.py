import os
import pygad
cwd = os.getcwd()
pygad_loc = '/grid/home/nbourgeois/codonOpt'
os.chdir(pygad_loc)
from general_functions import *
from metrics import *
os.chdir(cwd)
import time

        
def run_GA(aa_seq, 
           tissue, 
           filename, 
           population_size=20,
           cai_on = True, 
           bai_on = True, 
           cpg_on = True):
    
    # ftiness function
    def fitness_func(solution, solution_idx):

        if not type(solution) is str:
            seq_aa = ''.join([codon_to_int[x] for x in solution])
        else:
            seq_aa = solution
    #     print(solution_idx)

        tmp_dict = {}

        #Check for redundancy
        if seq_aa in all_sols.keys():
            fitness = all_sols[seq_aa]['fitness']

        else:
            fitness = 0

            if cai_on:
                cai = get_cai(seq_aa, cai_weight_dict)
                fitness += cai*cai_w
                tmp_dict['cai'] = cai

            if bai_on:
                bai = get_bai(seq_aa, bai_weight_dict)
                fitness += bai*bai_w
                tmp_dict['bai'] = bai

            if cpg_on:
                cpg = get_cpg(seq_aa)
                fitness += cpg*cpg_w
                tmp_dict['cpg'] = cpg

            fitness = fitness/total_weight
            tmp_dict['fitness'] = fitness
            all_sols[seq_aa] = tmp_dict

            
        return fitness

    def callback_generation(ga_instance):
        global last_fitness
#         print("Generation = {generation}".format(generation=ga_instance.generations_completed))
#         print("Fitness    = {fitness}".format(fitness=ga_instance.best_solution()[1]))
#         print("Change     = {change}".format(change=ga_instance.best_solution()[1] - last_fitness))
        last_fitness = ga_instance.best_solution()[1]
        toc = time.time()
        elapsed = toc-tic
        time_fit.append((elapsed,last_fitness))

    # Initialize
    time_fit = []
    tic = time.time()
    cai_w = 1 *cai_on
    bai_w = 1 *bai_on
    cpg_w = 1 *cpg_on
    total_weight = sum([cai_w,bai_w,cpg_w])
    cai_weight_dict = get_codon_weights(tissue)
    bai_weight_dict = get_bicodon_weights(tissue)
    codon_to_int, gene_space = init_parameters(aa_seq)
    gene_space_int = [[codon_to_int[x] for x in y] for y in gene_space]
    
    # fitness function variables
    desired_output = 1 # Function output.
    num_genes = len(gene_space)
    num_generations = 100 # Number of generations.
    num_parents_mating = 5 # Number of solutions to be selected as parents in the mating pool.
    sol_per_pop = population_size 
    parent_selection_type = "sss" # Type of parent selection.
    keep_parents = 5 # Number of parents to keep in the next population. -1 means keep all parents and 0 means keep nothing.
    crossover_type = "two_points" # Type of the crossover operator.
    # Parameters of the mutation operation.
    mutation_type = "random" # Type of the mutation operator.
    mutation_percent_genes = 5 # Percentage of genes to mutate. This parameter has no action if the parameter mutation_num_genes exists or when mutation_type is None.
    last_fitness = 0

    # Creating an instance of the GA class inside the ga module. Some parameters are initialized within the constructor.
    all_sols = {}
    fitness_function = fitness_func
    ga_instance = pygad.GA(num_generations=num_generations,
                           num_parents_mating=num_parents_mating, 
                           fitness_func=fitness_function,
                           sol_per_pop=sol_per_pop, 
                           num_genes=num_genes,
                           parent_selection_type=parent_selection_type,
                           keep_parents=keep_parents,
                           crossover_type=crossover_type,
                           mutation_type=mutation_type,
                           mutation_percent_genes=mutation_percent_genes,
                           on_generation=callback_generation,
                           gene_type=int,
                           gene_space=gene_space_int,
    )


    # Running the GA to optimize the parameters of the function.
    ga_instance.run()

    best_solution_generation=ga_instance.best_solution_generation

#     ga_instance.save(filename=filename)
    return(time_fit)

