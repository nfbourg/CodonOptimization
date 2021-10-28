import sys
import os 
from Bio import SeqIO
import pandas as pd

pygad_loc = os.path.dirname(os.path.abspath(__file__))

def tissue_opt(pygad_loc=pygad_loc):
    all_tiss = os.listdir(os.path.join(pygad_loc,'references/CoCoPUTs_codon_usage/codon_usage/'))
    tissues = [x.split('.')[0] for x in all_tiss]
    tissues.append('')
    tis_drop = widgets.Dropdown(
        options=tissues,
        value='',
        description='Tissue:',
        disabled=False,
    )
    return tis_drop

def readFasta( fastaFile: str ):
    fastaId = ''
    fastaSeq = ''
    if os.path.exists( fastaFile ):
        try:
            with open( fastaFile, "r") as handle:
                cont = 0
                fastaIds = []
                fastaSeqs = []
                for record in SeqIO.parse( handle, "fasta"):
                    fastaIds.append(record.id)
                    fastaSeqs.append(record.seq)
        except IOError as e:
            print( "I/O error({0}): {1}".format(e.errno, e.strerror) )
        except:
            print( "Unexpected error:", sys.exc_info()[0] )
    else:
        print("File does not exist:" + fastaFile )
        sys.exit(2)
    return (fastaIds, fastaSeqs)

def init_parameters(aa_seq, codon_usage_table_loc=os.path.join(pygad_loc,'references/codon_usage.getex.txt')):
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


    return(codon_to_int, gene_space)


