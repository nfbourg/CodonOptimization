import pandas as pd
import os
import glob
import mygene 
# os.chdir(os.path.join('tissue_cocoputs','data'))

# print(os.listdir()
# )
# human_cds = pd.read_csv('o586358-Human_CDS_Bicod.tsv')

# human_cds

def subset_and_write(df,col_list,fn):
    df = df.set_index('transcript_id')
    df = df[col_list]
    df.to_csv(fn)

def create_tissue_directory(sample_data, tissues):
    sample_dir = {}
    for tissue in tissues:
        sample_dir[tissue] = sample_data.loc[sample_data['SMTSD'] == tissue]['SAMPID'].values.tolist()
    return(sample_dir)

def get_weights(geneSyms, gtex_avg):

    def append_results(transripts,refseqs,weights):
        sub_transcripts = [x for x in transripts if x in gtex_avg.index]
        weight = gtex_avg[sub_transcripts].sum()
        refseqs.append(gdict['query'])
        weights.append(weight)
        
    refseqs = []
    weights = []
    counter = 0
    for gdict in geneSyms:
        if 'ensembl' in gdict:
            if len(gdict['ensembl']) == 1:
                transripts = gdict['ensembl']['transcript']
                append_results(transripts,refseqs,weights)
            else:
                # if the results has a list for the 'transcript' dict key
                for ensmbl in gdict['ensembl']:
                    transripts = ensmbl['transcript']
                    append_results(transripts,refseqs,weights)

        counter +=1    

    return(refseqs,weights)

def generate_bicodon_table(refseqs,weights,file_id):

    bcds_sorted = bcds.loc[refseqs,:]
    bcds_weighted = (bcds_sorted.T * weights).sum(axis=1)
    bcds_weighted.index = bcds_weighted.index.str.upper()
    bcds_weighted.to_frame().transpose().to_csv(
        'bicodon_tables/GTEx_{}_{}_v8.tsv'.format(tissue,file_id), sep='\t', index=False)
    print('Generated Bicodon Table')


def generate_codon_table(refseqs,weights,file_id):

    cds_sorted = cds.loc[refseqs,:]
    cds_weighted = (cds_sorted.T * weights).sum(axis=1)

    # output codon freq table
    output = pd.read_csv('/grid/home/nbourgeois/codonOpt/references/codon_usage.getex.txt',sep='\t')
    output = output.drop(columns='Fraction')
    output.loc[:,'Number'] = cds_weighted.values.astype(int)
    output.loc[:,'Frequency'] = round(output['Number'] / sum(output['Number']) * 1000,1)
    output.to_csv('codon_tables/GTEx_{}_{}_v8.tsv'.format(tissue,file_id), sep='\t', index=False)
    print('Generated Codon Table')

if __name__ == '__main__':

    ref_dir = os.path.join(os.getcwd(),'..','references')
    os.chdir('/grid/home/nbourgeois/codonOpt/tissue_cocoputs/data')

    # New Data
    gtex_data = pd.read_csv('ref/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct',sep='\t',skiprows=2)
    gtex_samp_data= pd.read_csv('ref/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt',sep='\t')

    print('Read in File')

    gtex_samp_data = gtex_samp_data.loc[gtex_samp_data['SAMPID'].isin(gtex_data.columns)]
        
    tissues = gtex_samp_data['SMTSD'].unique().tolist()
    sample_dir = create_tissue_directory(gtex_samp_data,tissues)

    # read in raw codon info
    cds = pd.read_csv('o586358-Human_CDS.tsv',sep='\t')
    cols = cds.columns
    cds = cds.drop(columns=['GGG'])
    cds = cds.reset_index()
    cds.columns = cols
    cds = cds.drop(columns = ['Gene ID', 'Accession', 'Division', 'Assembly', 'Taxid',
                            'Species', 'Organelle', '# Codons', 'GC%', 'GC1%', 'GC2%', 'GC3%'])
    cds = cds.groupby('Protein ID').sum()

    # read in raw codon info
    bcds = pd.read_csv('o586358-Human_CDS_Bicod.tsv',sep='\t')
    bcds = bcds.drop(columns = ['Gene ID', 'Accession', 'Division', 'Assembly', 'Taxid',
                            'Species', 'Organelle', '# Codon Pairs', 'GC%', 'GC1%', 'GC2%', 'GC3%','Unnamed: 4109'])
    bcds = bcds.groupby('Protein ID').sum()

    # convert protein info to ensmbl transcript info
    mg = mygene.MyGeneInfo()
    geneSyms = mg.querymany(cds.index.values, scopes='refseq',
                            fields='ensembl.transcript', species='human')


    for tissue in tissues:
        print(tissue)

        # subset to specific tissue
        cols = ['transcript_id','gene_id'] + sample_dir[tissue]
        gtex_data_tissue = gtex_data[cols]
        tissue = tissue.replace(" - ", "_").replace(" ", "_").replace("_(BA9)", "")

        # save tissue table
        gtex_data_tissue.to_csv('TPM/GTEx_{}_TPM_v8.tsv'.format(tissue),sep='\t')

        # reformat
        gtex_data_tissue['transcript_id'] = gtex_data_tissue['transcript_id'].str.split('.',expand=True)[0]
        gtex_data_tissue = gtex_data_tissue.groupby('transcript_id').sum()

        # get average per transcript
        gtex_avg = gtex_data_tissue.mean(axis=1)
        refseqs,weights = get_weights(geneSyms, gtex_avg)
        generate_codon_table(refseqs,weights,'codon_orig')
        generate_bicodon_table(refseqs,weights,'bicodon_orig')
        
        hifr_gtex_avg = gtex_avg[gtex_avg > 100]
        hifr_refseqs,hifr_weights = get_weights(geneSyms, hifr_gtex_avg)
        generate_codon_table(hifr_refseqs,hifr_weights,'codon_hifr')
        generate_bicodon_table(hifr_refseqs,hifr_weights,'bicodon_hifr')
        # get the liver weights for each protein id
        





