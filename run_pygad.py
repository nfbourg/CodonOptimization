import pygad
from general_functions import *
from metrics import *
import time
from multiprocessing.pool import Pool
from datetime import datetime


        
def run_GA(aa_seq, 
           tissue, 
           generations,
           cai_on = True, 
           bai_on = True, 
           cpg_on = True, 
           differential = True, 
           tissue2 = 'Heart_Atrial_Appendage',
           threads = 1):
    """Wrapper to initialize and run the genetic algorithm

    Args:
        aa_seq (str): amino acid string of the protein to be codon optimized
        tissue (str): tissue type of the protein for cai/bai.
        generations (int): Generations of the genetic algorithm
        cai_on (bool, optional): Turns cai metric on/off. Defaults to True.
        bai_on (bool, optional): Turns bai metric on/off. Defaults to True.
        cpg_on (bool, optional): Turns cpg metric on/off. Defaults to True.
        threads (int, optional): How many multiprocessing threads to use. Defaults to 1.

    
    Returns:
        GA_super: returns wrapped pygad.GA instance
    """
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")

    print(f'Start Time: {current_time}')
    pygad.GA = Updated_GA(tissue,cai_on,bai_on,cpg_on,differential,tissue2)

    ga_instance = GA_super(aa_seq, generations, threads)

    ga_instance.run()

    return ga_instance

class Updated_GA(pygad.GA):

    def __init__(self, tissue, cai_on, bai_on, cpg_on, differential=False,tissue2=False):
        """Add additional parameters to the default pygad class

        Args:
            pygad ([type]): [description]
            tissue ([type]): [description]
            cai_on ([type]): [description]
            cai_on ([type]): [description]
            cai_on ([type]): [description]

        Returns:
            [type]: [description]
        """

        pygad.GA.cai_on = cai_on
        pygad.GA.bai_on = bai_on
        pygad.GA.cpg_on = cpg_on
        pygad.GA.differential = differential
        pygad.GA.cai_weight_dict = get_codon_weights(tissue)
        pygad.GA.bai_weight_dict = get_bicodon_weights(tissue) 
        pygad.GA.bicodon_dict = {}

        if differential:
            pygad.GA.cai_weight_dict2 = get_codon_weights(tissue2)
            pygad.GA.bai_weight_dict2 = get_bicodon_weights(tissue2) 
            pygad.GA.bicodon_dict2 = {}


        pygad.GA.fit_dict = {}

class GA_super(pygad.GA):
    """Wrapper class for pygad that contains all additional initialization for the codon optimization"""
     
    def __init__(self, aa_seq, generations, threads):
        self.time_fit = []
        self.tic = time.time() # start timer
        self.gen = 1
        self.time_limit = 1e10 
        self.pool = Pool(processes=threads)    
        self.threads=threads
        self.init_parameters(aa_seq, generations)
        # set parameters on the pygad object 
        super().__init__(fitness_func=self.fitness_func,**self.variable_dict)

    def init_parameters(self, aa_seq, generations, 
        codon_usage_table_loc=os.path.join(pygad_loc,'references/codon_usage.getex.txt')):
        """Retrieve default parameters for pygad.  

        Pygad requires an input gene space that is integers, so the aa_seq needs to be converted to 
        an integer list with each integer corresponding to a unique codon. The function also requires
        allof the defualt parameters for pygad.

        Args:
            aa_seq (str): Input amino acid sequence to codon optimize
            codon_usage_table_loc (str): location of codon usage table
            
        Assigns:
            codon_to_int (dict): dict to map codons to ints and back
            gene_space_int (list[list]): list of lists where each index corresponds to a position
                in the input aa_seq, and the list in that index are all synonymous codons for the aa
                in the pos for the original seq
            variable_dict (dict): dictionary of following defualt values for pygad
                num_genes (int): number of genes in the gene space
                num_generations (int): Number of generations.
                num_parents_mating (int): Number of solutions to be selected as 
                    parents in the mating pool.
                sol_per_pop (int): population size 
                parent_selection_type (str): # Type of parent selection, 
                    sss is steady state selection 
                keep_parents (int): Number of parents to keep in the next population.
                        -1 means keep all parents and 0 means keep nothing.
                crossover_type (str): Type of the crossover operator. 
                mutation_type (str): Type of the mutation operator.
                mutation_percent_genes (int): Percentage of genes to mutate. 
                'stop_criteria': If the fitness fnction reaches 1, it will stop
                'gene_type' (type): expected object type for each gene.  Genes (In our 
                    case each codon) is uniquely mapped to an integer.
                'on_generation' (func): callback function to be called each gen
        """         

        codon_usage_table = pd.read_csv(codon_usage_table_loc,sep='\t')
        forward_table = pd.Series(codon_usage_table.AA.values,index=codon_usage_table.Codon).to_dict()

        back_table = {}
        for key in forward_table:
            val = forward_table[key]
            if val not in back_table.keys():
                back_table[val] = [key]
            else:
                back_table[val].append(key)

        codon_to_int = {}
        i=0
        for codon in forward_table.keys():
            codon_to_int[codon] = i
            codon_to_int[i] = codon
            i += 1

        gene_space = []
        for aa in aa_seq:
            all_cds = back_table[aa]
            gene_space.append(all_cds)

        gene_space_int = [[codon_to_int[x] for x in y] for y in gene_space]

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
                
        self.codon_to_int = codon_to_int
        self.gene_space_int = gene_space_int
        self.variable_dict = variables

    @staticmethod
    def fitness_func(solution, ind):
        """fitness function calculates metric for the 

        Args:
            solution ([list or str]): list of int or str of aa_seq
            ind (int): placeholder val to satisfy req of pygad class

        Returns:
            fitness (int): fitness value [0,1]
        """
        
        # weights for metrics
        cai_w = 1 *pygad.GA.cai_on
        bai_w = 1 *pygad.GA.bai_on
        cpg_w = 1 *pygad.GA.cpg_on
        total_weight = sum([cai_w,bai_w,cpg_w])

        fitness = 0

        # if not str, convert from list[int] to string
        if not type(solution) is str:
            codon_to_int = pygad.GA.codon_to_int
            seq_aa = ''.join([codon_to_int[x] for x in solution])
        else:
            seq_aa = solution

        if pygad.GA.cai_on:
            cai = get_cai(seq_aa, pygad.GA.cai_weight_dict)
            if pygad.GA.differential:
                cai2 = get_cai(seq_aa, pygad.GA.cai_weight_dict2)
                cai = cai-cai2
            fitness += cai*cai_w

        if pygad.GA.bai_on:
            # bai = get_bai(seq_aa, pygad.GA.bai_weight_dict, pygad.GA.bicodon_dict)
            bai = get_bai(seq_aa, pygad.GA.bai_weight_dict)
            if pygad.GA.differential:
                # bai2 = get_bai(seq_aa, pygad.GA.bai_weight_dict2, pygad.GA.bicodon_dict2)
                bai2 = get_bai(seq_aa, pygad.GA.bai_weight_dict2)
                bai = bai-bai2
            fitness += bai*bai_w

        if pygad.GA.cpg_on:
            cpg = get_cpg(seq_aa)
            fitness += cpg*cpg_w
            
        return fitness/total_weight
        
    def int_to_str(self,sol):
        """converts a codon list[int] to aa string"""
        seq_aa = ''.join([self.codon_to_int[x] for x in sol])
        return seq_aa

    def add_to_dict(self, seqs_to_add):
        """Add the sequences and their fitness to the fitness dict

        Args:
            seqs_to_add ([type]): [description]
        """
        if self.threads==1:
            pop_fitness = [self.fitness_func(x,0) for x in seqs_to_add]

        else:
            toc = time.time()
            elapsed = toc-self.tic
            if elapsed > self.time_limit:
                pop_fitness = [2 for x in seqs_to_add]
            else:
                data = [(x,0) for x in seqs_to_add]
                pop_fitness = self.pool.starmap(self.fitness_func, data)

        for i in range(len(seqs_to_add)): 
            self.fit_dict[seqs_to_add[i]] = pop_fitness[i]

    def cal_pop_fitness(self):
        """overwrites the base cal_pop_fitness function from pygad

        Rewrites the cal_pop_fitness to add a 1) a solution dictionary and 2) mp compatibility

        Returns:
            pop_fitness list[int]: list of population fitness
        """
        seqaa_population = [self.int_to_str(pop) for pop in self.population]
        seqs_missing = [x for x in seqaa_population if x not in self.fit_dict.keys()]
        self.add_to_dict(seqs_missing)
        pop_fitness = [self.fit_dict[x] for x in seqaa_population]
        pop_fitness = np.array(pop_fitness)
        return pop_fitness



def callback_generation(ga_instance):
    """callback function for each GA gen

    Args:
        ga_instance ([type]): [description]
    """
    if ga_instance.gen%100 == 0:
        print(f"Generation: {ga_instance.gen}")
        print(f"Time Passed: {time.time() - ga_instance.tic} seconds\n")
    ga_instance.gen = ga_instance.gen + 1

    last_fitness = ga_instance.best_solution()[1]
    # print("Generation = {generation}".format(generation=ga_instance.generations_completed))
    # print("Fitness    = {fitness}".format(fitness=ga_instance.best_solution()[1]))
    # print("Change     = {change}\n".format(change=ga_instance.best_solution()[1] - last_fitness))