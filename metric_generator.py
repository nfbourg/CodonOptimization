from general_functions import *
from metrics import *
from multiprocessing.pool import Pool
import pandas as pd
from collections import defaultdict
import RNA
from itertools import product
import numpy as np

class MetricGenerator():

    def __init__(self, tissue=None, fasta_loc=None, phe_data=None, wt_loc=None):

        if fasta_loc is not None:
            self.seq_ids,self.seqs = readFasta(fasta_loc)

        self.tissue=tissue
        self.func_dict = {}
        self.phe_data=phe_data
        if tissue is not None:
            self.init_cai_bai(tissue)
            self.init_rel_codon_freq(tissue)

        if phe_data is not None:
            self.phe_data = phe_data

        
        self.func_dict['CpG'] = get_cpg
        self.func_dict['Codon_Dist'] = get_codon_dist
        self.func_dict['Junc_Dist'] = count_cpg_junc
        self.func_dict['HSC'] = count_HSC
        self.func_dict['MFE'] = get_mfe
        self.func_dict['%CG'] = percent_cg
        self.func_dict['%CG3'] = percent_cg3
        self.init_junc_func()

        if wt_loc is not None:
            self.wt_seq = str(readFasta(wt_loc)[1][0])
            self.wt_metrics = {}
            for metric in self.func_dict.keys():
                self.wt_metrics[metric] = self.func_dict[metric](self.wt_seq)

    def init_rel_codon_freq(self,tissue):
        self.codon_freq = load_codon_freq_tab(tissue)
        self.func_dict['codon_freq'] = self.codon_freq_dist
        self.func_dict['codon_freq_sum'] = self.codon_freq_sum

    def init_cai_bai(self,tissue):
        self.cai_weight_dict = get_codon_weights(tissue)
        self.func_dict['CAI'] = self.get_cai
        self.func_dict['CAI-med'] = self.get_median_cai

        self.bai_weight_dict = get_bicodon_weights(tissue)
        self.func_dict['BAI'] = self.get_bai
        self.func_dict['BAI-percentile'] = self.get_percentile_bai

    def run_all_metrics(self):
        try: 
            return(self.all_metrics)
        except AttributeError as e:
            df = pd.DataFrame()

            for metric in self.func_dict.keys():
                df[metric] = self.run_metric(metric)
            df.index = self.seq_ids
            self.all_metrics = df
        if self.phe_data is None:
            return(df)
        else:
            try:
                df = df.merge(self.phe_data,left_index=True,right_index=True,how='left')
            except AttributeError as e:
                pass
            return(df)

    def run_metric(self, met_key):
        met_func = self.func_dict[met_key]
        metric_vals = [met_func(x) for x in self.seqs]
        return metric_vals

    def run_metric_arg(self, met_key, arg):
        met_func = self.func_dict[met_key]
        metric_vals = [met_func(x, arg) for x in self.seqs]
        return metric_vals

    def get_median_cai(self,seq):
        return get_median_cai(seq,self.cai_weight_dict)

    def get_cai(self,seq):
        return get_cai(seq,self.cai_weight_dict)

    def get_bai(self,seq):
        return get_bai(seq,self.bai_weight_dict)

    def get_percentile_bai(self,seq, percentile=10):
        return get_percentile_bai(seq,self.bai_weight_dict,percentile)

    def init_junc_func(self):
        nts = ['A','T','C','G']
        juncs = [x[0]+x[1] for x in product(nts, nts)]
        for junc in juncs:
            self.func_dict[f'Junc_freq_{junc}'] =  self.make_junc_func(junc)
        return
    
    def freq_range(self,seq):
        jdist = count_cpg_junc(seq)
        diff = jdist.max() - jdist.min()
        return diff

    def codon_freq_dist(self,seq):
        freqs = [self.codon_freq[seq[i:i+3]] for i in range(0, len(seq)-3, 3)] #convert to codon_list -> use weight dictionary
        return(freqs)

    def codon_freq_sum(self,seq):
        freqs = self.codon_freq_dist(seq)
        freqs = np.array(freqs) 
        freq_sum = freqs[freqs > .45].sum()
        return(freq_sum)
    # def calc_codon_freq():


    def make_junc_func(self,junc):

        def get_junc_freq(seq):
            jdist = count_cpg_junc(seq)
            if junc not in jdist.index:
                return(0)
            else:
                count = jdist[junc]
                freq = count/jdist.sum()
            return(freq)

        return get_junc_freq
    
    
def get_percentile_bai(seq,weight_dict,percentile):
    num_subseqs = int(round(100/percentile))
    bai_list = []

    if num_subseqs*3 > len(seq):
        print('Too small percentile')
        return([None]*num_subseqs)

    cds_per = int(np.ceil( (len(seq)/num_subseqs) / 3 ))
    start = 0
    while(start < len(seq)):
        stop = start + cds_per*3 + 3
        subseq = seq[start:stop]
        start = stop - 3
        bai = get_bai(subseq,weight_dict)
        bai_list.append(bai)

    return(bai_list)

def get_median_cai(seq,weight_dict):
    if type(weight_dict) is str:
        weight_dict = get_codon_weights(weight_dict)
    weights = [weight_dict[seq[i:i+3]] for i in range(0, len(seq)-3, 3)]
    med_cai = np.median(weights)
    return(med_cai)

def count_cpg_junc(seq):
    junctions = [i+j for i,j in zip(seq[2::3],seq[3::3])]
    j_dict = defaultdict(lambda: 0)
    for j in junctions:
        j_dict[j] += 1
    jdist=pd.Series(j_dict)
    return(jdist)

def count_HSC(seq):
    scs = ['TAG','TAA','TGA']
    counter = -1
    for sc in scs:
        counter += seq.count(sc)
    return(counter)

def calculate_MFE(seq, window):
    values=[]
    for i in range(len(seq)):
        temp_seq=seq[i:i+window]
        (ss, mfe) = RNA.fold(temp_seq)
        values.append(mfe)
    return(values)

def get_mfe(seq_aa):
    all_mfe = calculate_MFE(str(seq_aa),10)
    mfe = np.mean(all_mfe)  * -1
    return(mfe)

def percent_cg3(seq):
    seq = seq[2::3]
    cg = ['C','G']
    counter = 0
    for char in cg:
        counter += seq.count(char)
    return(counter/len(seq))

def percent_cg(seq):
    cg = ['C','G']
    counter = 0
    for char in cg:
        counter += seq.count(char)
    return(counter/len(seq))

def load_codon_freq_tab(tissue):
    # read table
    try:
        index_loc = os.path.join(pygad_loc,'references/CoCoPUTs_codon_usage/codon_usage/'+tissue+'.codon.txt')
        codon_table_raw = pd.read_csv(index_loc,sep='\t')

    except FileNotFoundError:
        print('Invalid Tissue')
        return()    

    codon_table = codon_table_raw.set_index('#Codon')
    sums = codon_table.groupby(by=['AA'])['Frequency'].sum()
    codon_table['sum_col'] = codon_table.apply(lambda x: sums[x['AA']],axis=1)
    codon_table['Rel_Freq'] = codon_table['Frequency'] / codon_table['sum_col']
        
    return(codon_table['Rel_Freq'])
# class OptSeq():


#     def __init__(self, tissue, nt_seq, generations):
#         """Add additional parameters to the default pygad class

#         Args:
#             tissue ([type]): [description]
#             aa_seq (str): Input amino acid sequence to codon optimize
#             codon_usage_table_loc (str): location of codon usage table
            
#         Assigns:
#             codon_to_int (dict): dict to map codons to ints and back
#             gene_space_int (list[list]): list of lists where each index corresponds to a position
#                 in the input aa_seq, and the list in that index are all synonymous codons for the aa
#                 in the pos for the original seq

#         Returns:
#             [type]: [description]
#         """

#         self.cai_weight_dict = get_codon_weights(tissue)
#         self.bai_weight_dict = get_bicodon_weights(tissue) 
#         self.fit_dict = {}
#         self.codon_to_int = codon_to_int
#         self.gene_space_int = gene_space_int
#         self.cai = get_cai()

#         codon_usage_table = pd.read_csv(codon_usage_table_loc,sep='\t')
#         forward_table = pd.Series(codon_usage_table.AA.values,index=codon_usage_table.Codon).to_dict()

#         back_table = {}
#         for key in forward_table:
#             val = forward_table[key]
#             if val not in back_table.keys():
#                 back_table[val] = [key]
#             else:
#                 back_table[val].append(key)

#         codon_to_int = {}
#         i=0
#         for codon in forward_table.keys():
#             codon_to_int[codon] = i
#             codon_to_int[i] = codon
#             i += 1

#         gene_space = []
#         for aa in aa_seq:
#             all_cds = back_table[aa]
#             gene_space.append(all_cds)

#         gene_space_int = [[codon_to_int[x] for x in y] for y in gene_space]
                

#         self.variable_dict = variables
        
#     def int_to_str(self,sol):
#         """converts a codon list[int] to aa string"""
#         seq_aa = ''.join([self.codon_to_int[x] for x in sol])
#         return seq_aa
