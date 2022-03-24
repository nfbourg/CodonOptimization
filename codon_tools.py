import os
import pandas as pd
import random 
from general_functions import *


pygad_loc = os.path.dirname(os.path.abspath(__file__))

def gen_random_seqs(prot_seq, num_seqs):
    codon_usage_table_loc = os.path.join(pygad_loc,'references','codon_usage.getex.txt')
    codon_to_int, codon_space = init_parameters_fun(prot_seq, codon_usage_table_loc) 
    
    positions = len(codon_space)
    size_vector = [len(codon_pos) for codon_pos in codon_space] 
    seqs = [None]*num_seqs
    for i in range(num_seqs):
        rand_index = [random.randrange(size) for size in size_vector]
        seqs[i] = ''.join(([space[ind] for ind, space in zip(rand_index,codon_space)]))
    return(seqs)

def translator(seq):
    global forward_table
    codon_usage_table_loc = os.path.join(pygad_loc,'references','codon_usage.getex.txt')
    codon_usage_table = pd.read_csv(codon_usage_table_loc,sep='\t')
    forward_table = pd.Series(codon_usage_table.AA.values,index=codon_usage_table.Codon).to_dict()
    
    back_table = {}
    for key in forward_table:
        val = forward_table[key]
        if val not in back_table.keys():
            back_table[val] = [key]
        else:
            back_table[val].append(key)
    
    for x,y in zip(seq_records[0],seq_records[1]):
        seq = str(y).upper()
        aa_seq = ''  
        for i in range(0, len(seq), 3):
            codon = seq[i:i+3]
            if len(codon) ==3:
                aa_seq = aa_seq + forward_table[codon]
            else:
                print('Sequence not divisible by 3, ends with:',codon)
    return(aa_seq)
