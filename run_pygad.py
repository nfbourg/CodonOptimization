import pygad
from general_functions import *
from metrics import *
import time
from multiprocessing.pool import Pool

        
def run_GA(aa_seq, 
           tissue, 
           generations,
           cai_on = True, 
           bai_on = True, 
           cpg_on = True,
           threads = 1):


    
    # Initialize
    cai_w = 1 *cai_on
    bai_w = 1 *bai_on
    cpg_w = 1 *cpg_on
    total_weight = sum([cai_w,bai_w,cpg_w])
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

    ga_instance = GA_super(codon_to_int,gene_space,gene_space_int,tissue,generations,threads)
    ga_instance = GA_super(codon_to_int,gene_space,gene_space_int,tissue,generations,threads)
    
    ga_instance.run()

    best_solution_generation=ga_instance.best_solution_generation

    return(ga_instance)

 
    


    # def fitness_wrapper(solution):
    #     return fitness_func(solution, 0)


class GA_super(pygad.GA):
    #     def __init__(self,codon_to_int,*args,**kwargs):
    #         self.pool = Pool(processes=2)    s
    #         self.all_sols = {}
    #         self.codon_to_int = codon_to_int
    # #         super().__init__(fitness_func=self.fitness_func,*args,**kwargs)
    #         variables = get_defaults()
    #         super().__init__(fitness_func=self.fitness_func,*args,**kwargs)


            
        def __init__(self,codon_to_int,gene_space,gene_space_int,tissue,generations,threads):
            self.time_fit = []
            self.tic = time.time()
            self.pool = Pool(processes=threads)    
            self.threads=threads
            pygad.GA.all_sols = {}
            pygad.GA.codon_to_int = codon_to_int
            pygad.GA.cai_weight_dict = get_codon_weights(tissue)
            pygad.GA.bai_weight_dict = get_bicodon_weights(tissue)            
            variable_dict = self.get_params(gene_space_int,generations)
            super().__init__(fitness_func=self.fitness_func,**variable_dict)

    
        # ftiness function
        @staticmethod
        def fitness_func(solution, ind):
            
            cai_w,bai_w,cpg_w = 1,1,1
            total_weight=3
            all_sols = {}
            if not type(solution) is str:
                codon_to_int = pygad.GA.codon_to_int
                seq_aa = ''.join([codon_to_int[x] for x in solution])
            else:
                seq_aa = solution

            tmp_dict = {}

            #Check for redundancy
            if seq_aa in all_sols.keys():
                fitness = all_sols[seq_aa]['fitness']

            else:
                fitness = 0

                cai_weight_dict = pygad.GA.cai_weight_dict
                bai_weight_dict = pygad.GA.bai_weight_dict
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
        def get_params(self,gene_space_int,generations):
            population_size = 750
            variables={
                    'num_generations':generations,
                    'sol_per_pop':population_size, 
                    'num_parents_mating':round(.25*population_size), 
                    'num_genes':len(gene_space_int),
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

        def int_to_str(self,sol):
            seq_aa = ''.join([self.codon_to_int[x] for x in sol])
            return seq_aa

        def add_to_dict(self, seqs_to_add):

            if self.threads==1:
                pop_fitness = [self.fitness_func(x,0) for x in seqs_to_add]

            else:
                toc = time.time()
                elapsed = toc-self.tic
                if elapsed > 10800:
                    pop_fitness = [2 for x in seqs_to_add]
                else:
                    data = [(x,0) for x in seqs_to_add]
                    pop_fitness = self.pool.starmap(self.fitness_func, data)

            for i in range(len(seqs_to_add)): 
                self.all_sols[seqs_to_add[i]] = pop_fitness[i]

        def cal_pop_fitness(self):
            seqaa_population = [self.int_to_str(pop) for pop in self.population]
            seqs_missing = [x for x in seqaa_population if x not in self.all_sols.keys()]
            self.add_to_dict(seqs_missing)
            pop_fitness = [self.all_sols[x] for x in seqaa_population]
            pop_fitness = np.array(pop_fitness)
            return pop_fitness

def callback_generation(ga_instance):
#         print("Generation = {generation}".format(generation=ga_instance.generations_completed))
#         print("Fitness    = {fitness}".format(fitness=ga_instance.best_solution()[1]))
#         print("Change     = {change}".format(change=ga_instance.best_solution()[1] - last_fitness))
    # global stop
    last_fitness = ga_instance.best_solution()[1]
    # toc = time.time()
    # elapsed = toc-ga_instance.tic
    # ga_instance.time_fit.append((elapsed,last_fitness))