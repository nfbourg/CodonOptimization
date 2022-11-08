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

    def __init__(self, aa_seq, tissues, ntissues=None, negative=False,
                 prefix_codon=None, ramp=False, wt_ramp=False, wt_seq='', mimic=False, cpg_max=None):
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
        self.mimic = mimic
        self.negative = negative
        self.ramp = ramp
        self.wt_ramp = wt_ramp
        self.wt_seq = wt_seq
        self.cpg_max = cpg_max

        if ramp or wt_ramp:
            self.init_ramp(wt_seq)

        if mimic:
            self.init_mimic(depth=1)

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

    def init_mimic(self,depth):
        self.depth=depth

    def init_ramp(self,wt_seq):

        self.start_bai = .4
        self.min_bai  = .2
        self.ramp_end = 600
        self.ramp_start = 350

        if self.wt_ramp:
            aa_start = int( self.ramp_start/3)
            for i in range(aa_start):
                ind = i*3
                self.gene_space[i] = [wt_seq[ind:ind+3]]

    def calculate_tissue_bai(self,seq):
        bais = []
        for tissue in self.tissues:
            bais.append(get_bai(seq,self.bai_weight_dict[tissue]))
        return geo_mean(bais)    


    # What is this for???
    # def trim(self, nt_ind):
    #     codon_ind = int(np.floor((nt_ind-1)/3))
    #     codon = self.result[codon_ind]
    #     if len(self.gene_space)>1:
    #         self.gene_space[codon_ind].remove(codon)
    #         for i in range(codon_ind+1, len(self.result)):
    #             self.codon_chains[i] = None
    #     else:
    #         print('Cannot trim, gene space too small.')
    #     return self.optimize()

    def optimize(self):
        self.check_chain(len(self.result)-1)
        self.result = self.codon_chains[-1]['TAA'].copy()

        if self.prefix is not None:
            end_result = self.result[1:]
        else:
            end_result = self.result

        if self.added:
            final_seq = ''.join(end_result[:-1])
            print(self.calculate_tissue_bai(final_seq))
            return(final_seq)
        else:
            final_seq = ''.join(end_result)
            print(self.calculate_tissue_bai(final_seq))
            return(final_seq)


    def check_chain(self,ind):

        if self.codon_chains[ind] is None:
            self.check_chain(ind-1)
        print(ind,end=' ')
        self.chain_codon(ind)

    def score_chain(self, chain):

        def calc_chain_bai(chain, tissue):
            self.perf_count+=1
            seq = ''.join(chain) + 'TAA'
            if self.method == 'CAI':
                sub_bai = get_cai(seq,self.cai_weight_dict[tissue])
            else:
                sub_bai = get_bai(seq,self.bai_weight_dict[tissue])
            return(sub_bai)

        if not self.mimic:
            sub_bai_list = [calc_chain_bai(chain,tissue) for tissue in self.tissues]
            sub_bai_gmean = geo_mean(sub_bai_list)
        

        if self.differential:
            sub_neg_bai_list = [calc_chain_bai(chain,tissue) for tissue in self.ntissues]
            sub_bai_gmean = sub_bai_gmean - geo_mean(sub_neg_bai_list)
                        ## somehow zeroes appear below


        elif self.mimic:

            depth = self.depth + 2 # default dyn program needs to go back 2 for bai 

            new_bai_list = [calc_chain_bai(chain[-depth:],tissue) for tissue in self.tissues]

            new_bai_mean = geo_mean(new_bai_list)    


            wt_start = max(0,len(chain)*3 - depth*3)
            wt_end = len(chain)*3 
            wt_bai_list = [get_bai(self.wt_seq[wt_start:wt_end]+ 'TAA',self.bai_weight_dict[tissue]) for tissue in self.tissues]
            wt_bai_gmean = geo_mean(wt_bai_list)    

            bai_target = wt_bai_gmean + .2
            bai_target = min(bai_target,1)

            sub_bai_gmean = 1 - abs(new_bai_mean - bai_target)
            # print(new_bai_mean, wt_bai_gmean)

            if new_bai_mean < wt_bai_gmean: # discourage dropping below wt
                sub_bai_gmean = .01
                # return(sub_bai_gmean*.5)

        elif self.ramp:
            if len(chain)*3 <= self.ramp_end: 

                sub_bai_list = [calc_chain_bai(chain[-20:],tissue) for tissue in self.tissues]
                sub_bai_gmean = geo_mean(sub_bai_list)    

                if sub_bai_gmean < self.min_bai:
                    return(0)

                ramp_target = 1 - (1-self.start_bai) * (self.ramp_end -len(chain)*3 +1)/self.ramp_end 
                sub_bai_gmean = 1 - abs(sub_bai_gmean-ramp_target)

        elif self.wt_ramp:
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


        if self.negative:
            return(1-sub_bai_gmean)

        else:
            return(sub_bai_gmean)

    def cpg_adjustment(self, chain, chain_score):

        def calc_chain_cpg(chain):

            depth = self.depth + 2

            seq = ''.join(chain[-depth:])
            sub_cpg = 1-get_cpg(seq)

            return(sub_cpg)
       
        # check for cpgs
        if self.cpg_max is not None:
            cpg_perc = calc_chain_cpg(chain)
            if cpg_perc <= self.cpg_max:
                pass
            else:
                divisor = (cpg_perc + .001) / ( + .001)
                chain_score = chain_score * cpg_perc / divisor
        
        return(chain_score)

    def chain_codon(self, ind):

        self.codon_chains[ind] = {}
        for codon in self.gene_space[ind]:
            cmax=-999

            if ind+1 < len(self.gene_space): # cehck if last codon
                peek_ind = ind+1
                max_ind = min(peek_ind+self.depth, len(self.gene_space))
                new_codon_list = product([codon],*self.gene_space[peek_ind:max_ind])
            else: 
                new_codon_list = [[codon]]

            for new_codons in new_codon_list:
                for chain_key in self.codon_chains[ind-1]:
                    new_chain = self.codon_chains[ind-1][chain_key].copy()
                    new_chain.extend(new_codons)
                    chain_score = self.score_chain(new_chain)
                    chain_score = self.cpg_adjustment(new_chain, chain_score)

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

