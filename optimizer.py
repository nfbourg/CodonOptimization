from general_functions import *
from itertools import product
from metrics import *


class Optimizer:

    def __init__(
        self, aa_seq, opt_type, tissues, negative=False, method='BAI', **kwargs
        ):
        """_summary_

        Args:
            aa_seq (_type_): _description_
            tissues (_type_): _description_
            negative (bool, optional): True if negatively optimizing. Defaults to False.
            ntissues (_type_, optional): If negatively opti. Defaults to None.
            opt_type (str, optional): Must be in [max, target, ramp, wt_ramp, mimic]. Defaults to max.
            prefix_codon (_type_, optional): _description_. Defaults to None.
            wt_seq (str, optional): _description_. Defaults to ''.
            cpg_thresh (_type_, optional): _description_. Defaults to None.

        Optional Args:
            negative
            cpg_thresh
            wt_seq
            prefix_codon
            ntissues
        """
        def get_value(key):
            try:
                wt_seq = kwargs[key]
            except KeyError:
                raise(f'Missing parameter: {key}')
            return(wt_seq)

        def input_to_list(tissues):
            if type(tissues) is str:
                return [tissues]
            else:
                return list(tissues)

        self.tissues = input_to_list(tissues)#, add * to beginning of input arg

        if not aa_seq.endswith('*'):
            aa_seq = aa_seq + '*'
            self.added=True
        else:
            self.added=False

    
        self.method = method
        self.negative = negative

        try:
            self.prefix_codon = kwargs['prefix_codon']
            self.result = [None] * ( len(aa_seq)+1 )
            self.result[0] = self.prefix_codon

        except KeyError:
            self.prefix_codon = None
            self.result = [None] * len(aa_seq)
            self.result[0] = 'ATG'

        try:
            self.ramp_end = int(kwargs['ramp_end'])
        except KeyError:
            self.ramp_end = 600

            
        try:
            self.ramp_start = int(kwargs['ramp_start'])
        except KeyError:
            self.ramp_start = 350


        self.init_tissues(self.tissues)

        self.init_parameters(aa_seq)
        self.init_codon_chains()
        self.opt_type = opt_type
        match opt_type:

            case 'max':
                # required args: [method, ]
                pass

            case 'differential':
                self.ntissues = input_to_list(get_value('ntissues'))
                self.init_tissues(self.ntissues, init=False)

                # required args: [method, ]

            case 'target':
                # required args: [method, target]
                self.target = get_value('target')

            case 'ramp':
                # required args: [method, ]
                self.wt_ramp = False
                self.init_ramp(depth=3)


            case 'wt_ramp':
                # required args: [method, ]
                self.wt_seq = get_value('wt_seq')
                self.wt_ramp = True
                self.init_ramp(depth=3)

            case 'mimic':
                # required args: [method, ]
                self.wt_seq = get_value('wt_seq')
                self.target_range = get_value('target_range')
                self.depth = get_value('depth')
                # self.init_mimic()

            case _:
                raise('Incorrect optimization type')



        if 'cpg_thresh' in kwargs.keys():
            self.cpg_thresh=kwargs['cpg_thresh']

        else:
            self.cpg_thresh=None



    def init_codon_chains(self):
        self.codon_chains = [None] * len(self.result)
        self.codon_chains[0] = {self.result[0]:[self.result[0]]}
        self.codon_chains[1] = {}
        for codon in self.gene_space[1]:
            self.codon_chains[1][codon] = [self.result[0],codon]



    def init_cpg(self,cpg_thresh):
        self.cpg_thresh = cpg_thresh
        self.cpg_min = True

    def init_tissues(self,tissues,init=True):

        if init:
            self.bai_weight_dict = {}
            self.cai_weight_dict = {}

        # self.tissues = input_to_list(tissues)
        # dict_tissues = self.tissues.copy()

        # if ntissues is not None:
        #     self.ntissues = input_to_list(ntissues)
        #     dict_tissues.extend(self.ntissues)

        for tissue in tissues:
            self.bai_weight_dict[tissue] = get_bicodon_weights(tissue) 
            self.cai_weight_dict[tissue] = get_codon_weights(tissue) 

    def init_mimic(self,depth=1,target_range=.4):
        self.depth=depth
        self.target_range=target_range


    def init_ramp(self,depth):

        self.start_bai = .4
        self.min_bai  = .2
        self.depth=depth

        if self.wt_ramp:
            aa_start = int( self.ramp_start/3)
            for i in range(aa_start):
                ind = i*3
                self.gene_space[i] = [self.wt_seq[ind:ind+3]]

    def calculate_tissue_bai(self,seq):
        bais = []
        for tissue in self.tissues:
            bais.append(get_bai(seq,self.bai_weight_dict[tissue]))
        return geo_mean(bais)    


    def optimize(self):
        self.check_chain(len(self.result)-1)
        self.result = self.codon_chains[-1]['TAA'].copy()

        if self.prefix_codon is not None:
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

        match self.opt_type:

            case 'max' | 'target' | 'differential':

                sub_bai_list = [calc_chain_bai(chain,tissue) for tissue in self.tissues]
                sub_bai_gmean = geo_mean(sub_bai_list)
           
                match self.opt_type:

                    case 'differential':
                        sub_neg_bai_list = [calc_chain_bai(chain,tissue) for tissue in self.ntissues]
                        sub_bai_gmean = sub_bai_gmean - geo_mean(sub_neg_bai_list)
                                    ## somehow zeroes appear below
                    case 'target':
                        sub_bai_gmean = 1 - np.sqrt(abs(((sub_bai_gmean)**2 - (self.target)**2)))
                    
                    case 'max':
                        pass


            case 'mimic':

                depth = self.depth + 2 # default dyn program needs to go back 2 for bai 

                new_bai_list = [calc_chain_bai(chain[-depth:],tissue) for tissue in self.tissues]

                new_bai_mean = geo_mean(new_bai_list)    

                wt_start = max(0,len(chain)*3 - depth*3)
                wt_end = len(chain)*3 
                wt_bai_list = [get_bai(self.wt_seq[wt_start:wt_end]+ 'TAA',self.bai_weight_dict[tissue]) for tissue in self.tissues]
                wt_bai_gmean = geo_mean(wt_bai_list)    

                bai_target = wt_bai_gmean + self.target_range
                bai_target = min(bai_target,1)

                sub_bai_gmean = 1 - np.sqrt(abs(((new_bai_mean)**2 - (bai_target)**2)))
                # print(new_bai_mean, wt_bai_gmean)

                if new_bai_mean < wt_bai_gmean: # discourage dropping below wt
                    # sub_bai_gmean = .0001
                    sub_bai_gmean = (1- (wt_bai_gmean-new_bai_mean))*.01

            case 'ramp':
                if len(chain)*3 <= self.ramp_end: 

                    depth = self.depth + 2
                    sub_bai_list = [calc_chain_bai(chain[-depth:],tissue) for tissue in self.tissues]
                    sub_bai_gmean = geo_mean(sub_bai_list)    

                    if sub_bai_gmean < self.min_bai:
                        return(0)

                    ramp_target = 1 - (1-self.start_bai) * (self.ramp_end -len(chain)*3 +1)/self.ramp_end 
                    sub_bai_gmean = 1 - abs(sub_bai_gmean-ramp_target)
                else:
                    sub_bai_list = [calc_chain_bai(chain,tissue) for tissue in self.tissues]
                    sub_bai_gmean = geo_mean(sub_bai_list)

            case 'wt_ramp':
                if len(chain)*3 < self.ramp_start:
                    sub_bai_list = [calc_chain_bai(chain,tissue) for tissue in self.tissues]
                    sub_bai_gmean = geo_mean(sub_bai_list)
                    self.start_bai = sub_bai_gmean
                    return(1)
                if len(chain)*3 > self.ramp_end: 
                    sub_bai_list = [calc_chain_bai(chain,tissue) for tissue in self.tissues]
                    sub_bai_gmean = geo_mean(sub_bai_list)
                else:
                    depth = self.depth + 2

                    sub_bai_list = [calc_chain_bai(chain[-depth:],tissue) for tissue in self.tissues]
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
        if self.cpg_thresh is not None:
            cpg_perc = calc_chain_cpg(chain)
            if cpg_perc <= self.cpg_thresh:
                pass
            else:
                divisor = (cpg_perc + .001) / ( + .001)
                chain_score = chain_score * cpg_perc / divisor
        
        return(chain_score)

    def chain_codon(self, ind):
        """Optimizing at position ind, with the assumption that position ind-1 
        is already has the most optimal chain calculated for all codons in ind-1 
        codon space of ind

        Args:
            ind (int): index that will be optimized
        """
        def remove_cpg(new_codon_list):
            new_codon_list = [codons for codons in new_codon_list if 'CG' not in ''.join(codons)]
            return(new_codon_list)

        self.codon_chains[ind] = {}
        for codon in self.gene_space[ind]:
            cmax=-999
            
            if ind+1 < len(self.gene_space): # cehck if last codon
                peek_ind = ind+1
                max_ind = min(peek_ind+self.depth, len(self.gene_space))
                new_codon_list = product([codon],*self.gene_space[peek_ind:max_ind])
            else: 
                new_codon_list = [[codon]]

            # if self.no_cpg:
            #     new_codon_list = remove_cpg(new_codon_list)

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


    def init_parameters(self, aa_seq, 
            codon_usage_table_loc=os.path.join(pygad_loc,'references/codon_usage.getex.txt')):
        
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

        if self.prefix_codon is not None:
            gene_space.insert(0,self.prefix_codon)

        # gene_space_int = [[codon_to_int[x] for x in y] for y in gene_space]
                
        self.gene_space = gene_space

        self.perf_count=0
        self.depth=1
        # self.codon_to_int = codon_to_int
        # self.gene_space_int = gene_space_int

    def optimizer_main(self):

     
        optimized_seq = self.optimize()
        seq_name = f"mimic_wt_{target_range}d{depth}"
        


        match self.opt_type:
            
            case 'mimic':
                depth = params['depth']
                target_range = params['target_range']
                self.init_mimic(depth,target_range)
        
            case 'ramp':
                depth = params['depth']
                seq_name = f"ramp_d{depth}"
                self.depth=depth
                optimized_seq = self.optimize()

            case 'wt_ramp':
                depth = params['depth']
                self.depth=depth
                optimized_seq = self.optimize()
                seq_name = f"wt_ramp_d{depth}"

            case 'cpg':
                depth = params['depth']
                thresh = params['cpg_thresh']
                tar=.85
                self.target=tar
                self.depth=depth
                optimized_seq = self.optimize()
                
                seq_name = f"cpg{round(cpg_perc,2)}_tar{tar}d{depth}t{thresh}"
            
            case _ :
                print('ERROR')
                
        fileout = f"vectors/{seq_name}.fa"
        with open(fileout, "w") as fileo:
            fileo.write(f">{seq_name}\n{optimized_seq}\n")

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