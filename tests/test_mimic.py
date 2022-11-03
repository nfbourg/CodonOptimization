

import os
cwd = os.getcwd()
import sys
sys.path.insert(1, '/grid/home/nbourgeois/codonOpt')
from general_functions import *
import dyn_prog
import importlib


wt_loc = '/grid/home/nbourgeois/data/codon_jason/wt_pah.fa'
a,seq = readFasta(wt_loc)
wt_seq = str(seq[0])
ga_input = '/grid/home/nbourgeois/data/test_proteins/pah/pah.pep.fas' #sequence pep.fas Input
tissue = 'Liver' # Tissue type for CoCoPuts

(keys, seqs) = readFasta(ga_input)
if len(seqs) == 1:
    aa_seq=str(seqs[0])

aa_seq = 'MSTAVLENP'
wt_seq  = wt_seq[:len(aa_seq)*3]

importlib.reload(dyn_prog)
optimizer = dyn_prog.Optimizer(aa_seq, tissues='Liver', mimic=True, wt_seq=wt_seq)
optimized_seq = optimizer.optimize()
