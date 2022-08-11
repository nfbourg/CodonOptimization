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

    def __init__(self, aa_seq, tissues='Liver', ntissues=None):
        if not aa_seq.endswith('*'):
            aa_seq = aa_seq + '*'
            self.added=True
        else:
            self.added=False

        if type(tissues) is str:
            self.tissues = [tissues]
        else:
            self.tissues = list(tissues)

        self.bai_weight_dict = {}
        for tissue in self.tissues:
            self.bai_weight_dict[tissue] = get_bicodon_weights(tissue) 

        if ntissues is None:
            self.differential=False
        else:
            if type(ntissues) is str:
                self.ntissues = [ntissues]
            else:
                self.ntissues = list(ntissues)

            for ntissue in self.ntissues:
                self.bai_weight_dict[ntissue] = get_bicodon_weights(ntissue) 
                self.differential=True

        self.init_parameters(aa_seq)
        self.depth=1

        self.result = [None] * len(aa_seq)
        self.result[0] = 'ATG'
        self.codon_chains = [None] * len(aa_seq)
        self.codon_chains[0] = {self.result[0]:[self.result[0]]}
        self.codon_chains[1] = {}
        for codon in self.gene_space[1]:
            self.codon_chains[1][codon] = [self.result[0],codon]

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
        if self.added:
            return(''.join(self.result[:-1]))
        else:
            return(''.join(self.result))

   
    def check_chain(self,ind):
        # for codon1 in self.gene_space[ind-2]:
        if self.codon_chains[ind] is None:
            self.check_chain(ind-1)

        self.chain_codon(ind)
        print(ind,end=' ')



    def chain_codon(self, ind):

        def calc_chain_bai(new_chain,tissue):
            seq = ''.join(new_chain) + 'TAA'
            sub_bai = get_bai(seq,self.bai_weight_dict[tissue])
            return(sub_bai)

        self.codon_chains[ind] = {}
        for codon in self.gene_space[ind]:
            cmax=-999
            if ind+1 < len(self.gene_space):
                peek_ind = ind+1
                max_ind = min(peek_ind+self.depth, len(self.gene_space))
                new_codon_products = product([codon],*self.gene_space[peek_ind:max_ind])
                for new_codons in new_codon_products:
                # for codon2 in self.gene_space[ind+1]:
                    # new_codons = [codon,codon2]
                    for chain_key in self.codon_chains[ind-1]:
                        # print(new_codons)
                        new_chain = self.codon_chains[ind-1][chain_key].copy()
                        new_chain.extend(new_codons)
                                    
                        sub_bai_list = [calc_chain_bai(new_chain,tissue) for tissue in self.tissues]
                        sub_bai_gmean = geo_mean(sub_bai_list)

                        if self.differential:
                            sub_neg_bai_list = [calc_chain_bai(new_chain,tissue) for tissue in self.ntissues]
                            sub_bai_gmean = sub_bai_gmean - geo_mean(sub_neg_bai_list)
                        # if (0 in sub_bai_list):                 ## somehow zeroes appear below
                        #     sub_bai_gmean = -2
                        # else:
                        if sub_bai_gmean > cmax:
                            cmax = sub_bai_gmean
                            self.codon_chains[ind][codon] = self.codon_chains[ind-1][chain_key].copy()
                            self.codon_chains[ind][codon].append(codon)
            else:
                for chain_key in self.codon_chains[ind-1]:
                    new_chain = self.codon_chains[ind-1][chain_key].copy()
                    new_chain.append(codon)
                                
                    sub_bai_list = [calc_chain_bai(new_chain,tissue) for tissue in self.tissues]
                    # if (0 in sub_bai_list):                 ## somehow zeroes appear below
                    #     sub_bai_gmean = -2
                    # else:
                    sub_bai_gmean = geo_mean(sub_bai_list)
                    if sub_bai_gmean > cmax:
                        cmax = sub_bai_gmean
                        self.codon_chains[ind][codon] = new_chain
                                            # if self.differential:
                    #     tri_neg_bai_list = [calc_tri_bai(bicodon1,bicodon2,ntissue) for ntissue in self.ntissues]
                    #     if (0 in tri_bai_list):                 ## somehow zeroes appear below
                    #         tri_bai_gmean = -2
                    #     else:
                    #         tri_bai_gmean = tri_bai_gmean - geo_mean(tri_neg_bai_list)
        
        # return(maxcodon1,maxcodon2,maxcodon3)


    # def get_maximal_bai(self, ind):

    #     def calc_tri_bai(bicodon1,bicodon2,tissue):
    #         tri_bai = np.exp(len( self.result) * np.log( self.bai_weight_dict[tissue][bicodon1] * self.bai_weight_dict[tissue][bicodon2] ) )
    #         return(tri_bai)

    #     cmax=-999
    #     codon1 = self.result[ind-2]
    #     for codon2 in self.gene_space[ind-1]:
    #         for codon3 in self.gene_space[ind]:
    #             bicodon1 = codon1 + codon2
    #             bicodon2 = codon2 + codon3

    #             tri_bai_list = [calc_tri_bai(bicodon1,bicodon2,tissue) for tissue in self.tissues]
    #             if (0 in tri_bai_list):                 ## somehow zeroes appear below
    #                 tri_bai_gmean = -2
    #             else:
    #                 tri_bai_gmean = geo_mean(tri_bai_list)

    #             if self.differential:
    #                 tri_neg_bai_list = [calc_tri_bai(bicodon1,bicodon2,ntissue) for ntissue in self.ntissues]
    #                 if (0 in tri_bai_list):                 ## somehow zeroes appear below
    #                     tri_bai_gmean = -2
    #                 else:
    #                     tri_bai_gmean = tri_bai_gmean - geo_mean(tri_neg_bai_list)

    #             if tri_bai_gmean > cmax:
    #                 cmax=tri_bai_gmean   
    #                 maxcodon1 = codon1
    #                 maxcodon2 = codon2
    #                 maxcodon3 = codon3
    #     return(maxcodon1,maxcodon2,maxcodon3)

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

        # gene_space_int = [[codon_to_int[x] for x in y] for y in gene_space]
                
        self.gene_space = gene_space
        # self.codon_to_int = codon_to_int
        # self.gene_space_int = gene_space_int

