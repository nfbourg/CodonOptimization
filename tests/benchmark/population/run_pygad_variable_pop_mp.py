import os
import pygad
cwd = os.getcwd()
pygad_loc = '/grid/home/nbourgeois/codonOpt'
os.chdir(pygad_loc)
from general_functions import *
from metrics import *
os.chdir(cwd)
import time
from multiprocessing.pool import Pool
import functools  

        
def run_GA(aa_seq, 
           tissue, 
           population_size=20,
           cai_on = True, 
           bai_on = True, 
           cpg_on = True,
           threads = 1):
    
    # Initialize
    # time_fit = []
    # tic = time.time()
    global stop 
    stop = False
    cai_w = 1 *cai_on
    bai_w = 1 *bai_on
    cpg_w = 1 *cpg_on
    total_weight = sum([cai_w,bai_w,cpg_w])
#     cai_weight_dict = get_codon_weights(tissue)
#     bai_weight_dict = get_bicodon_weights(tissue)
    codon_to_int, gene_space = init_parameters(aa_seq)
    gene_space_int = [[codon_to_int[x] for x in y] for y in gene_space]
    
#     ========================
#     # fitness function variables
#     desired_output = 1 # Function output.
#     num_genes = len(gene_space)
#     num_generations = 100 # Number of generations.
#     num_parents_mating = 5 # Number of solutions to be selected as parents in the mating pool.
#     sol_per_pop = population_size 
#     parent_selection_type = "sss" # Type of parent selection.
#     keep_parents = 5 # Number of parents to keep in the next population. -1 means keep all parents and 0 means keep nothing.
#     crossover_type = "two_points" # Type of the crossover operator.
#     # Parameters of the mutation operation.
#     mutation_type = "random" # Type of the mutation operator.
#     mutation_percent_genes = 5 # Percentage of genes to mutate. This parameter has no action if the parameter mutation_num_genes exists or when mutation_type is None.
#     ========================

    last_fitness = 0


#     all_sols = {}
#     fitness_function = fitness_func
                
#     ga_instance = GA_super(codon_to_int=codon_to_int,
#                            num_generations=num_generations,
#                            num_parents_mating=num_parents_mating, 
# #                            fitness_func=fitness_function,
#                            sol_per_pop=sol_per_pop, 
#                            num_genes=num_genes,
#                            parent_selection_type=parent_selection_type,
#                            keep_parents=keep_parents,
#                            crossover_type=crossover_type,
#                            mutation_type=mutation_type,
#                            mutation_percent_genes=mutation_percent_genes,
#                            on_generation=callback_generation,
#                            gene_type=int,
#                            gene_space=gene_space_int,
#     )

#     with Pool(processes=5) as pool:
        # Running the GA to optimize the parameters of the funcabstion.
    ga_instance = GA_super(codon_to_int,gene_space,gene_space_int,population_size,tissue,threads)
    
    ga_instance.run()


    best_solution_generation=ga_instance.best_solution_generation

    return(ga_instance.time_fit)

def wraps_partial(wrapper, *args, **kwargs):
    """ Creates a calslable object whose attributes will be set from the partials nested func attribute ..."""
    wrapper = wrapper.func
    while isinstance(wrapper, functools.partial):
        wrapper = wrapper.func
    return functools.wraps(wrapper, *args, **kwargs)

class GA_super(pygad.GA):
#     def __init__(self,codon_to_int,*args,**kwargs):
#         self.pool = Pool(processes=2)    s
#         self.all_sols = {}
#         self.codon_to_int = codon_to_int
# #         super().__init__(fitness_func=self.fitness_func,*args,**kwargs)
#         variables = get_defaults()
#         super().__init__(fitness_func=self.fitness_func,*args,**kwargs)


        
    def __init__(self,codon_to_int,gene_space,gene_space_int,population_size,tissue,threads):
        self.time_fit = []
        self.tic = time.time()
        self.pool = Pool(processes=threads)    
        self.threads=threads
        pygad.GA.all_sols = {}
        pygad.GA.codon_to_int = codon_to_int
        pygad.GA.cai_weight_dict = get_codon_weights(tissue)
        pygad.GA.bai_weight_dict = get_bicodon_weights(tissue)
        # self.fitpart=wraps_partial(functools.partial(self.fitness_func,
        #                 codon_to_int=codon_to_int,
        #                 cai_weight_dict=self.cai_weight_dict,
        #                 bai_weight_dict=self.bai_weight_dict
        #                ))
#         self.fitpart = lambda x,y: self.fitness_func(x,y, 
#                         codon_to_int=codon_to_int,
#                         cai_weight_dict=self.cai_weight_dict,
#                         bai_weight_dict=self.bai_weight_dict)
        
        variable_dict = self.get_defaults(gene_space,gene_space_int,population_size)
        super().__init__(fitness_func=self.fitness_func,**variable_dict)

    # ftiness function
    @staticmethod
    def fitness_func(solution, ind):
        codon_to_int = pygad.GA.codon_to_int
        cai_weight_dict = pygad.GA.cai_weight_dict
        bai_weight_dict = pygad.GA.bai_weight_dict
        
        cai_w,bai_w,cpg_w = 1,1,1
        total_weight=3
        all_sols = {}
        if not type(solution) is str:
            seq_aa = ''.join([codon_to_int[x] for x in solution])
        else:
            seq_aa = solution

        tmp_dict = {}

        #Check for redundancy
        if seq_aa in all_sols.keys():
            fitness = all_sols[seq_aa]['fitness']

        else:
            fitness = 0

#             if cai_on:
            cai = get_cai(seq_aa, cai_weight_dict)
            fitness += cai*cai_w
            tmp_dict['cai'] = cai

#             if bai_on:
            bai = get_bai(seq_aa, bai_weight_dict)
            fitness += bai*bai_w
            tmp_dict['bai'] = bai

#             if cpg_on:
            cpg = get_cpg(seq_aa)
            fitness += cpg*cpg_w
            tmp_dict['cpg'] = cpg

            fitness = fitness/total_weight
            tmp_dict['fitness'] = fitness
            all_sols[seq_aa] = tmp_dict

            
        return fitness
        
    # @staticmethod
    def get_defaults(self,gene_space,gene_space_int,population_size):
        variables={
                'num_generations':1000000,
                'num_parents_mating':round(.25*population_size), 
                'sol_per_pop':population_size, 
                'num_genes':len(gene_space),
                'parent_selection_type':"sss",
                'keep_parents':round(.25*population_size),
                'crossover_type':"two_points",
                'mutation_type':"random",
                'mutation_percent_genes':5,
                'on_generation':callback_generation,
                'gene_type':int,
                'stop_criteria':['reach_1'],
                'gene_space':gene_space_int}
        return(variables)

    def cal_pop_fitness(self):
        if self.threads==1:
            pop_fitness = [self.fitness_func(x,0) for x in self.population]

        else:
            toc = time.time()
            elapsed = toc-self.tic
            if elapsed > 600:
                pop_fitness = [2 for x in self.population]
            else:
                data = [(x,0) for x in self.population]
                pop_fitness = self.pool.starmap(self.fitness_func, data)
            # self.all_sols[x] = y for x,y in zip(self.population,pop_fitness)
        
        pop_fitness = np.array(pop_fitness)
        return pop_fitness
    


    # def fitness_wrapper(solution):
    #     return fitness_func(solution, 0)


def callback_generation(ga_instance):
#         print("Generation = {generation}".format(generation=ga_instance.generations_completed))
#         print("Fitness    = {fitness}".format(fitness=ga_instance.best_solution()[1]))
#         print("Change     = {change}".format(change=ga_instance.best_solution()[1] - last_fitness))
    global stop
    last_fitness = ga_instance.best_solution()[1]
    toc = time.time()
    elapsed = toc-ga_instance.tic
    ga_instance.time_fit.append((elapsed,last_fitness))