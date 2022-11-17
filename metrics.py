import gc
import time
import getopt
import sys
# import hashlib
import os
import subprocess as sp
import collections
import pandas as pd
from pkg_resources import resource_filename
# from spliceai.utils import one_hot_encode
import numpy as np

pygad_loc = os.path.dirname(os.path.abspath(__file__))

def get_cpg(seq):
    count = seq.count("CG")
    cpg_score = ( len(seq) - count * 2 ) / len(seq)  
    return(cpg_score)

# def predictSites_multi_broken( fastaIds,fastaSeqs,threshold=.1):

#     gpu0_avail = gpu_memory_usage(0)
#     gpu1_avail = gpu_memory_usage(1)
#     if gpu0_avail + gpu1_avail > 30000:
#         print('Out of memory\n')
#         return(-1)
#     if gpu0_avail < gpu1_avail:
#         gpu="0"
#     else:
#         gpu="1"
    
#     os.environ["CUDA_VISIBLE_DEVICES"]=gpu

#     from keras.models import load_model
#     import tensorflow as tf
#     from tensorflow import keras

#     results = {}
#     gpus =tf.config.experimental.list_physical_devices('GPU')

#     try: 
#         tf.config.experimental.set_virtual_device_configuration(gpus[0],
#             [tf.config.experimental.VirtualDeviceConfiguration(memory_limit=5000)])
        
#         context = 10000
#         keras.backend.set_learning_phase(0)
#         paths = ('models/spliceai{}.h5'.format(x) for x in range(1, 6))
#         models = [load_model(resource_filename('spliceai', x), compile=False ) for x in paths]
#         for seq_id, seq in zip(fastaIds,fastaSeqs):

#             x = one_hot_encode('N'*(context//2) + str(seq) + 'N'*(context//2))[None, :]
#             y = np.mean([models[m].predict(x) for m in range(5)], axis=0)
#             acceptor_prob = y[0, :, 1]
#             donor_prob = y[0, :, 2]
#             coord = list( range(1, len( donor_prob ) ) )
                    
#             # custom metric
#             results[seq_id] = (len(acceptor_prob[acceptor_prob > threshold]) + len(donor_prob[donor_prob > threshold]))/2
#         return(results)
#     except tf.errors as e:
#         # Splice AI custome sequeence prediction
#         print(e,'\n')
#         return(-1)

def predictSites_multi( fastaIds,fastaSeqs,threshold=.1):

    gpu0_avail = gpu_memory_usage(0)
    gpu1_avail = gpu_memory_usage(1)
    if gpu0_avail + gpu1_avail > 30000:
        print('Out of memory\n')
        return(-1)
    if gpu0_avail < gpu1_avail:
        gpu="0"
    else:
        gpu="1"
    
    os.environ["CUDA_VISIBLE_DEVICES"]=gpu

    from keras.models import load_model
    import tensorflow as tf
    from tensorflow import keras

    results = []
    gpus =tf.config.experimental.list_physical_devices('GPU')

    try: 
        tf.config.experimental.set_virtual_device_configuration(gpus[0],
            [tf.config.experimental.VirtualDeviceConfiguration(memory_limit=5000)])
        
        context = 10000
        keras.backend.set_learning_phase(0)
        paths = ('models/spliceai{}.h5'.format(x) for x in range(1, 6))
        models = [load_model(resource_filename('spliceai', x), compile=False ) for x in paths]
        for seq_id, seq in zip(fastaIds,fastaSeqs):

            x = one_hot_encode('N'*(context//2) + str(seq) + 'N'*(context//2))[None, :]
            y = np.mean([models[m].predict(x) for m in range(5)], axis=0)
            acceptor_prob = y[0, :, 1]
            donor_prob = y[0, :, 2]
            coord = list( range(1, len( donor_prob ) ) )
                    
            # custom metric
#             results[seq_id] = (len(acceptor_prob[acceptor_prob > threshold]) + len(donor_prob[donor_prob > threshold]))/2
            acceptor_count = 0
            donor_count = 0
            for i in coord:
                if donor_prob[i] > threshold :
                    donor_count += donor_prob[i]
                if acceptor_prob[i] > threshold and donor_count > 0:
                    acceptor_count += acceptor_prob[i]
            # custom metric
            results[seq_id] = acceptor_count * donor_count
        return(results)
    except tf.errors as e:
        # Splice AI custome sequeence prediction
        print(e,'\n')
        return(-1)

def gpu_memory_usage(gpu_id):
    command = f"nvidia-smi --id={gpu_id} --query-gpu=memory.used --format=csv"
    output_cmd = sp.check_output(command.split())
    memory_used = output_cmd.decode("ascii").split("\n")[1]
    # Get only the memory part as the result comes as '10 MiB'
    memory_used = int(memory_used.split()[0])

    return memory_used

def get_codon_weights(tissue):
    # read table
    try:
        index_loc = os.path.join(pygad_loc,'references/CoCoPUTs_codon_usage/codon_usage/'+tissue+'.codon.txt')
        codon_table_raw = pd.read_csv(index_loc,sep='\t')

    except FileNotFoundError:
        print('Invalid Tissue')
        return()    

    codon_table = codon_table_raw.set_index('#Codon')
    maxes = codon_table_raw.groupby(by=['AA'])['Frequency'].max()
    codon_table['max_col'] = codon_table.apply(lambda x: maxes[x['AA']],axis=1)
    codon_table['weight'] = codon_table['Frequency'] / codon_table['max_col']
    weight_dict = codon_table['weight'].to_dict()
    
    return(weight_dict)

def get_cai(seq, weight_dict):
    if type(weight_dict) is str:
        weight_dict = get_codon_weights(weight_dict)
    weights = [weight_dict[seq[i:i+3]] for i in range(0, len(seq)-3, 3)] #convert to codon_list -> use weight dictionary
    cai = geo_mean(weights)
    return(cai)


def get_codon_dist(seq):
    'calculate the distribution of how many of each codons there are'
    codons = [str(seq[i:i+3]) for i in range(0,len(seq),3) ]
    codon_dict = collections.defaultdict(lambda: 0)
    for codon in codons:
        codon_dict[codon] += 1/len(seq)
    return(codon_dict)


def geo_mean(iterable):
    a = np.array(iterable)
    a = np.log(a)
    return np.exp(a.sum() / len(a))

def get_bicodon_weights(tissue):
    
    def bicodon_to_AA(string):
        cod1 = string[:3]
        cod2 = string[3:]
        biAA = forward_table[cod1] + forward_table[cod2]
        return(biAA)
    
    def get_forward_table():
        codon_usage_table_loc = os.path.join(pygad_loc,'references/codon_usage.getex.txt')
        codon_usage_table = pd.read_csv(codon_usage_table_loc,sep='\t')
        forward_table = pd.Series(codon_usage_table.AA.values,index=codon_usage_table.Codon).to_dict()
        return(forward_table)
     
    # read table
    try:
        index_loc = os.path.join(pygad_loc,'references/CoCoPUTs_codon_usage/bicodon_usage/',tissue+'.bicodon.txt')
        codon_table_raw_t = pd.read_csv(index_loc,sep='\t')

    except FileNotFoundError:
        print('Invalid Tissue')
        return()
    
    codon_table = codon_table_raw_t.T

    # add amino acid pair column
    bicodons_nt = list(codon_table.index)
    forward_table = get_forward_table()
    codon_table['AApair'] = list(map(bicodon_to_AA,bicodons_nt))
    total = codon_table[0].sum()
    
    #add frequency
    bicodons_nt = list(codon_table.index)
    codon_table['AApair'] = list(map(bicodon_to_AA,bicodons_nt))
    total = codon_table[0].sum()
    codon_table['Frequency'] = codon_table.apply(lambda x: (x[0] / total * 1000) ,axis=1)
    
    #calculate weights
    maxes = codon_table.groupby(by=['AApair'])['Frequency'].max()
    codon_table['max_col'] = codon_table.apply(lambda x: maxes[x['AApair']],axis=1)
    codon_table['weight'] = codon_table['Frequency'] / codon_table['max_col']
    
    weight_dict = codon_table['weight'].to_dict()
    
    return(weight_dict)

def get_bai(seq, weight_dict):
    if type(weight_dict) is str:
        weight_dict = get_bicodon_weights(weight_dict)
       
    
    # sliding window bicodon, excludes the last codon
    weights = [weight_dict[seq[i:i+6]] for i in range(0, len(seq)-6, 3)] #convert to codon_list -> use weight dictionary
    
    bai = geo_mean(weights)
    return(bai)

def get_bai_sped(seq, weight_dict, bicodon_dict):
    if type(weight_dict) is str:
        weight_dict = get_bicodon_weights(weight_dict)
        
    slicer_size = 10
    ind=0
    weights = []
    seq = seq[:-3]
    while ind < (len(seq)):
        sub_seq = seq[ind:ind+(slicer_size*3)]
        if sub_seq in bicodon_dict.keys():
            tmp_weight = bicodon_dict[sub_seq]
        else:
            tmp_weights = [weight_dict[sub_seq[i:i+6]] for i in range(0, len(sub_seq)-3, 3)]
            print(tmp_weights)
            tmp_weight = np.prod(tmp_weights)
            bicodon_dict[sub_seq] = tmp_weight
        weights.append(tmp_weight)
        ind+=(slicer_size-1)*3
    # sliding window bicodon, excludes the last codon
#     weights = [weight_dict[seq[i:i+6]] for i in range(0, len(seq)-6, 3)] #convert to codon_list -> use weight dictionary
    tmp = np.array(weights)
    tmp = np.log(tmp)
    bai = np.exp(tmp.sum() / (len(seq)/3 -1))
#     bai = (np.prod(weights))**(3/(len(seq)))
    return(bai)
# def fitness_func(solution, solution_idx):
    
#     global all_sols
    
#     if not type(solution) is str:
#         seq_aa = ''.join([codon_to_int[x] for x in solution])
#     else:
#         seq_aa = solution
# #     print(solution_idx)

#     tmp_dict = {}
    
#     #Check for redundancy
#     if seq_aa in all_sols.keys():
#         fitness = all_sols[seq_aa]['fitness']

#     else:
#         fitness = 0
        
#         if cai_on:
#             cai = get_cai(seq_aa, cai_weight_dict)
#             fitness += cai*cai_w
#             tmp_dict['cai'] = cai
        
#         if bai_on:
#             bai = get_bai(seq_aa, bai_weight_dict)
#             fitness += bai*bai_w
#             tmp_dict['bai'] = bai
            
#         if cpg_on:
#             cpg = get_cpg(seq_aa)
#             fitness += cpg*cpg_w
#             tmp_dict['cpg'] = cpg

#         if sps_on:
#             sps = get_sps(seq_aa)
#             print('SPS retuned.')

#             fitness += sps*sps_w
#             tmp_dict['sps'] = sps

#         if pas_on:
#             pas = get_pas(seq_aa)
#             fitness += pas*pas_w
#             tmp_dict['pas'] = pas

#         fitness = fitness/total_weight
#         tmp_dict['fitness'] = fitness
#         all_sols[seq_aa] = tmp_dict
        
    
#     return fitness

