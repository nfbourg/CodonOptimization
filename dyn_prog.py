from concurrent.futures import thread
from re import L
from general_functions import *
from metrics import *
import time
from multiprocessing.pool import Pool
from datetime import datetime

class Optimizer():

    def __init__(self, aa_seq):

        self.bai_weight_dict = get_bicodon_weights('Liver') 
        self.init_parameters(aa_seq)
        self.result = list(range(len(aa_seq)))
        self.standby = {}


    def optimize(self,ind):
        if ind==0:
            self.result[ind] = 'ATG'
            self.codons = 'ATG'
            return(1)
            # for codon in self.codon_space[ind]:
            #     bicodon = 'ATG' + codon
            #     bai = (self.bai_weight_dict['ATG' + codon])
            #     if bai > cmax:
            #         cmax=bai
            #         self.result[ind] = codon
            return(cmax)
        else:
            self.optimize(ind-1)

            for codon in self.codon_space[ind]:
                bicodon = self.result[ind-1] + codon
                bai = (self.bai_weight_dict[bicodon])
                if bai > cmax:
                    cmax=bai
                    self.result[ind] = codon
            



            last = self.optimize(ind-1)
            lcodon = self.result[ind-1]

    def init_parameters(self, aa_seq,
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
                    
            self.gene_space = gene_space
            self.codon_to_int = codon_to_int
            self.gene_space_int = gene_space_int

