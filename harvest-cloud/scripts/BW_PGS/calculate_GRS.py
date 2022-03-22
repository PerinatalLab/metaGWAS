import pandas as pd
import numpy as np

def remove_palindromic(d, REF, EFF):
        return(d.loc[~(((d[REF]== 'T') & (d[EFF]== 'A')) | ((d[REF]== 'A') & (d[EFF]== 'T')) | ((d[REF]== 'C') & (d[EFF]== 'G')) | ((d[REF]== 'G') & (d[EFF]== 'C'))), :])

def flip_alleles(x):
        x= x.str.upper()
        x= x.str.replace('C', 'g')
        x= x.str.replace('G', 'c')
        x= x.str.replace('T', 'a')
        x= x.str.replace('A', 't')
        return(x.str.upper())

def harmonize_alleles(d, REF_x, EFF_x, REF_y, EFF_y):
#        d= remove_palindromic(d, REF_x, EFF_x)
#        d= remove_palindromic(d, REF_y, EFF_y)
        d['beta']= np.where((d[EFF_y]== d[EFF_x]) & (d[REF_y] == d[REF_x]), d.beta, np.where((d[EFF_y]== d[REF_x]) & (d[REF_y]== d[EFF_x]), -1 * d.beta, np.where((flip_alleles(d[EFF_y])== d[EFF_x]) & (flip_alleles(d[REF_y])== d[REF_x]), d.beta, np.where((flip_alleles(d[EFF_y])== d[REF_x]) & (flip_alleles(d[REF_y]) == d[EFF_x]), -1 * d.beta, np.nan))))
        d= d.loc[~(d.beta.isnull()), :]
        d[EFF_y]= np.where((d[EFF_y]== d[EFF_x]) & (d[REF_y]== d[REF_x]), d[EFF_x], np.where((d[EFF_y]== d[REF_x]) & (d[REF_y] == d[EFF_x]), d[REF_x], np.where((flip_alleles(d[EFF_y])== d[EFF_x]) & (flip_alleles(d[REF_y])== d[REF_x]), d[EFF_x], np.where((flip_alleles(d[EFF_y])== d[REF_x]) & (flip_alleles(d[REF_y]) == d[EFF_x]), d[REF_x], np.nan))))
        d[REF_y]= np.where((d[EFF_y]== d[EFF_x]) & (d[REF_y]== d[REF_x]), d[REF_x], np.where((d[EFF_y]== d[REF_x]) & (d[REF_y]== d[EFF_x]), d[EFF_x], np.where((flip_alleles(d[EFF_y])== d[EFF_x]) & (flip_alleles(d[REF_y])== d[REF_x]), d[REF_x], np.where((flip_alleles(d[EFF_y])== d[REF_x]) & (flip_alleles(d[REF_y]) == d[EFF_x]), d[EFF_x], np.nan))))
        return(d)

def calculate_PGS(d, ID):
        d.drop_duplicates(['chr', 'pos'], keep=False, inplace= True)
        d['chr']= d.chr.apply(str)
        d= pd.merge(betas, d, on= ['chr', 'pos'])
        d= harmonize_alleles(d, 'ref', 'eff', 'REF', 'EFF')
        d.drop_duplicates(['chr', 'pos'], keep=False, inplace= True)
        betas_v= np.array(d.beta)
        d.drop(['ref', 'eff', 'pos', 'chr', 'EFF', 'REF', 'beta'], axis= 1, inplace= True)
        ids= d.columns
        d= pd.DataFrame(np.array(d) * betas_v.reshape(-1, 1), columns= ids)
        d= pd.DataFrame(d.sum(axis= 0))
        d[ID]= d.index
        return d


# Read data
if 'haplotype' not in snakemake.input[0]:
        cols= ['chr','pos','ref','eff'] + [line.strip() for line in open(snakemake.input[0], 'r')]
        betas= pd.read_csv(snakemake.input[2], sep= '\t', header= 0)
        betas['chr']= betas['chr'].apply(str)
        df_list= list()
        for d in pd.read_csv(snakemake.input[1], header= None, names= cols, sep= '\t', chunksize= 600):
                d= calculate_PGS(d, 'IID')
                df_list.append(d)
        d= pd.concat(df_list)
else:
        betas= pd.read_csv(snakemake.input[1], sep= '\t', header= 0)
        betas['chr']= betas['chr'].apply(str)
        df_list= list()
        for d in pd.read_csv(snakemake.input[0], header= 0, sep= '\t', chunksize= 600):
                d= calculate_PGS(d, 'PREG_ID')
                df_list.append(d)
        d= pd.concat(df_list)

if 'haplotype' in snakemake.input[0]:
        d= d.groupby('PREG_ID').sum().reset_index()
        d.columns= ['PREG_ID', snakemake.wildcards.haplo]
else:
        d= d.groupby('IID').sum().reset_index()
        d.columns.values[0]= snakemake.wildcards.sample + '_' + snakemake.wildcards.repr_trait + '_PGS'


d.to_csv(snakemake.output[0], sep ='\t', header= True, index= False)

