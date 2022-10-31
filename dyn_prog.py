# from concurrent.futures import thread
# from re import L
from collections import defaultdict
from general_functions import *
from itertools import product
from metrics import *
# import time
# from multiprocessing.pool import Pool
# from datetime import datetime

class Optimizer():

    def __init__(self, aa_seq, tissues, ntissues=None, 
                 prefix_codon=None, ramp=False, og_ramp=False, og_seq=''):
        if not aa_seq.endswith('*'):
            aa_seq = aa_seq + '*'
            self.added=True
        else:
            self.added=False

        self.init_tissues(tissues,ntissues)
    
        self.init_parameters(aa_seq,prefix_codon)

        self.prefix = prefix_codon

        # ramp params
        self.method = 'BAI'
        self.ramp = ramp
        self.og_ramp = og_ramp
        if ramp or og_ramp:
            self.init_ramp(og_seq)

        if prefix_codon is None:
            self.result = [None] * len(aa_seq)
            self.result[0] = 'ATG'
        else:
            self.result = [None] * ( len(aa_seq)+1 )
            self.result[0] = prefix_codon

        self.codon_chains = [None] * len(self.result)
        self.codon_chains[0] = {self.result[0]:[self.result[0]]}
        self.codon_chains[1] = {}
        for codon in self.gene_space[1]:
            self.codon_chains[1][codon] = [self.result[0],codon]

    def init_tissues(self,tissues,ntissues):

        def input_to_list(tissues):
            if type(tissues) is str:
                return [tissues]
            else:
                return list(tissues)

        self.tissues = input_to_list(tissues)
        dict_tissues = self.tissues.copy()

        if ntissues is None:
            self.differential=False
        else:
            self.differential=True
            self.ntissues = input_to_list(ntissues)
            dict_tissues.extend(self.ntissues)

        self.bai_weight_dict = {}
        self.cai_weight_dict = {}
        for tissue in dict_tissues:
            self.bai_weight_dict[tissue] = get_bicodon_weights(tissue) 
            self.cai_weight_dict[tissue] = get_codon_weights(tissue) 

    def init_ramp(self,og_seq):
        self.start_bai = .4
        self.min_bai  = .2
        self.ramp_end = 600
        self.ramp_start = 350

        if self.og_ramp:
            aa_start = int( self.ramp_start/3)
            for i in range(aa_start):
                ind = i*3
                self.gene_space[i] = [og_seq[ind:ind+3]]

    def calculate_tissue_bai(self,seq):
        bais = []
        for tissue in self.tissues:
            bais.append(get_bai(seq,self.bai_weight_dict[tissue]))
        return geo_mean(bais)    

    def trim(self, nt_ind):
        codon_ind = int(np.floor((nt_ind-1)/3))
        codon = self.result[codon_ind]
        if len(self.gene_space)>1:
            self.gene_space[codon_ind].remove(codon)
            for i in range(codon_ind+1, len(self.result)):
                self.codon_chains[i] = None
        else:
            print('Cannot trim, gene space too small.')
        return self.optimize()

    def optimize(self):
        self.check_chain(len(self.result)-1)
        self.result = self.codon_chains[-1]['TAA'].copy()
        if self.prefix is not None:
            end_result = self.result[1:]
        else:
            end_result = self.result
        if self.added:
            return(''.join(end_result[:-1]))
        else:
            return(''.join(end_result))

    def check_chain(self,ind):
        # for codon1 in self.gene_space[ind-2]:
        if self.codon_chains[ind] is None:
            self.check_chain(ind-1)

        self.chain_codon(ind)

    def score_chain(self, chain):

        def calc_chain_bai(chain, tissue):
            self.perf_count+=1
            seq = ''.join(chain) + 'TAA'
            if self.method == 'CAI':
                sub_bai = \
            else:
                sub_bai = get_bai(seq,self.bai_weight_dict[tissue])
            return(sub_bai)
  

        sub_bai_list = [calc_chain_bai(chain,tissue) for tissue in self.tissues]
        sub_bai_gmean = geo_mean(sub_bai_list)

        if self.differential:
            sub_neg_bai_list = [calc_chain_bai(chain,tissue) for tissue in self.ntissues]
            sub_bai_gmean = sub_bai_gmean - geo_mean(sub_neg_bai_list)
                        ## somehow zeroes appear below
                        #     sub_bai_gmean = -2
                        # else:

        if self.ramp:
            if len(chain)*3 > self.ramp_end: 
                return(sub_bai_gmean)
            else:

                sub_bai_list = [calc_chain_bai(chain[-20:],tissue) for tissue in self.tissues]
                sub_bai_gmean = geo_mean(sub_bai_list)    

                if sub_bai_gmean < self.min_bai:
                    return(0)

                ramp_target = 1 - (1-self.start_bai) * (self.ramp_end -len(chain)*3 +1)/self.ramp_end 
                sub_bai_gmean = 1 - abs(sub_bai_gmean-ramp_target)

        if self.og_ramp:
            if len(chain)*3 < self.ramp_start:
                self.start_bai = sub_bai_gmean
                return(1)
            if len(chain)*3 > self.ramp_end: 
                return(sub_bai_gmean)
            else:

                sub_bai_list = [calc_chain_bai(chain[-10:],tissue) for tissue in self.tissues]
                sub_bai_gmean = geo_mean(sub_bai_list)    
                
                if sub_bai_gmean < self.min_bai:
                    return(0)

                ramp_target = 1 - (1-self.start_bai) * (self.ramp_end -len(chain)*3 +1)/(self.ramp_end - self.ramp_start)
                sub_bai_gmean = 1 - abs(sub_bai_gmean-ramp_target)
        


        

        return(sub_bai_gmean)

    def chain_codon(self, ind):

        self.codon_chains[ind] = {}
        for codon in self.gene_space[ind]:
            cmax=-999

            if ind+1 < len(self.gene_space): # cehck if last codon
                peek_ind = ind+1
                max_ind = min(peek_ind+self.depth, len(self.gene_space))
                new_codon_list = product([codon],*self.gene_space[peek_ind:max_ind])
            else: 
                new_codon_list = [codon]

            for new_codons in new_codon_list:
                for chain_key in self.codon_chains[ind-1]:
                    new_chain = self.codon_chains[ind-1][chain_key].copy()
                    new_chain.extend(new_codons)
                    chain_score = self.score_chain(new_chain)

                    if chain_score > cmax:
                        cmax = chain_score
                        self.codon_chains[ind][codon] = self.codon_chains[ind-1][chain_key].copy()
                        self.codon_chains[ind][codon].append(codon)

    def init_parameters(self, aa_seq, prefix_codon=None,
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

        # codon_to_int = {}
        # i=0
        # for codon in forward_table.keys():
        #     codon_to_int[codon] = i
        #     codon_to_int[i] = codon
        #     i += 1

        gene_space = []
        for aa in aa_seq:
            all_cds = back_table[aa]
            gene_space.append(all_cds)

        if prefix_codon is not None:
            gene_space.insert(0,prefix_codon)

        # gene_space_int = [[codon_to_int[x] for x in y] for y in gene_space]
                
        self.gene_space = gene_space

        self.perf_count=0
        self.depth=1
        # self.codon_to_int = codon_to_int
        # self.gene_space_int = gene_space_int

