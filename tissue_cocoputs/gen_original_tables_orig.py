import pandas as pd
import os
import glob
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

if __name__ == '__main__':

    ref_dir = os.path.join(os.getcwd(),'..','references')
    os.chdir('/grid/home/nbourgeois/codonOpt/tissue_cocoputs/data')

    # New Data
    # gtex_data = pd.read_csv('GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct',sep='\t',skiprows=2)
    # gtex_samp_data= pd.read_csv('GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt',sep='\t')

    # Original Data
    gtex_data = pd.read_csv('ref/GTEx_Analysis_2016-01-15_v7_RSEMv1.2.22_transcript_tpm.txt',sep='\t')
    gtex_samp_data= pd.read_csv('ref/GTEx_v7_Annotations_SampleAttributesDS.txt',sep='\t')

    print('Read in File')# gtex_data = pd.read_csv('GTEx_head',sep='\t',skiprows=2)


    gtex_samp_data = gtex_samp_data.loc[gtex_samp_data['SAMPID'].isin(gtex_data.columns)]
        
    tissues = gtex_samp_data['SMTSD'].unique().tolist()
    sample_dir = create_tissue_directory(gtex_samp_data,tissues)
    tissues = ['Liver']

    for tissue in tissues:
        print(tissue)
        # cols = ['gene_id'] + sample_dir[tissue]
        cols = ['transcript_id','gene_id'] + sample_dir[tissue]
        gtex_data_tissue = gtex_data[cols]
        tissue = tissue.replace(" - ", "_").replace(" ", "_").replace("_(BA9)", "")
        # gtex_data_tissue.to_csv('GTEx_{}_TPM.tsv'.format(tissue),sep='\t')
        gtex_data_tissue.to_csv('GTEx_{}_TPM_orig.tsv'.format(tissue),sep='\t')

        # gtex_data_tissue['gene_id'] = gtex_data_tissue['gene_id'].str.split('.',expand=True)[0]
        # gtex_data_tissue = gtex_data_tissue.groupby('gene_id').sum()
        # tissue_weight = gtex_data_tissue.mean(axis=1)
        # tissue_weight.to_csv('GTEx_{}_weight.tsv'.format(tissue),sep='\t')
        # print('Done')

   



