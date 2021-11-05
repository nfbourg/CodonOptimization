from general_functions import *
from metrics import *

def fitness_func(solution, solution_idx):
    
    global all_sols
    
    if not type(solution) is str:
        seq_aa = sol_to_str(solution)
    else:
        seq_aa = solution
#     print(solution_idx)

    tmp_dict = {}
    
    #Check for redundancy
    if seq_aa in all_sols.keys():
        fitness = all_sols[seq_aa]['fitness']

    else:
        fitness = 0
        
        if cai_on:
            cai = get_cai(seq_aa, cai_weight_dict)
            fitness += cai*cai_w
            tmp_dict['cai'] = cai
        
        if bai_on:
            bai = get_bai(seq_aa, bai_weight_dict)
            fitness += bai*bai_w
            tmp_dict['bai'] = bai
            
        if cpg_on:
            cpg = get_cpg(seq_aa)
            fitness += cpg*cpg_w
            tmp_dict['cpg'] = cpg

        fitness = fitness/total_weight
        tmp_dict['fitness'] = fitness
        all_sols[seq_aa] = tmp_dict
        
    
    return fitness

    
def fitness_func_orig(solution, solution_idx):
    
    global all_sols
    
    if not type(solution) is str:
        seq_aa = ''.join([codon_to_int[x] for x in solution])
    else:
        seq_aa = solution
#     print(solution_idx)

    tmp_dict = {}
    
    #Check for redundancy
    if seq_aa in all_sols.keys():
        fitness = all_sols[seq_aa]['fitness']

    else:
        fitness = 0
        
        if cai_on:
            cai = get_cai(seq_aa, cai_weight_dict)
            fitness += cai*cai_w
            tmp_dict['cai'] = cai
        
        if bai_on:
            bai = get_bai(seq_aa, bai_weight_dict)
            fitness += bai*bai_w
            tmp_dict['bai'] = bai
            
        if cpg_on:
            cpg = get_cpg(seq_aa)
            fitness += cpg*cpg_w
            tmp_dict['cpg'] = cpg

        if sps_on:
            sps = get_sps(seq_aa)
            print('SPS retuned.')

            fitness += sps*sps_w
            tmp_dict['sps'] = sps

        if pas_on:
            pas = get_pas(seq_aa)
            fitness += pas*pas_w
            tmp_dict['pas'] = pas

        fitness = fitness/total_weight
        tmp_dict['fitness'] = fitness
        all_sols[seq_aa] = tmp_dict
        
    
    return fitness

    