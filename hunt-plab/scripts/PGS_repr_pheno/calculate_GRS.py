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
        d= remove_palindromic(d, REF_x, EFF_x)
        d= remove_palindromic(d, REF_y, EFF_y)
        d['beta']= np.where((d[EFF_y]== d[EFF_x]) & (d[REF_y] == d[REF_x]), d.beta, np.where((d[EFF_y]== d[REF_x]) & (d[REF_y]== d[EFF_x]), -1 * d.beta, np.where((flip_alleles(d[EFF_y])== d[EFF_x]) & (flip_alleles(d[REF_y])== d[REF_x]), d.beta, np.where((flip_alleles(d[EFF_y])== d[REF_x]) & (flip_alleles(d[REF_y]) == d[EFF_x]), -1 * d.beta, np.nan))))
        d= d.loc[~(d.beta.isnull()), :]
        d[EFF_y]= np.where((d[EFF_y]== d[EFF_x]) & (d[REF_y]== d[REF_x]), d[EFF_x], np.where((d[EFF_y]== d[REF_x]) & (d[REF_y] == d[EFF_x]), d[REF_x], np.where((flip_alleles(d[EFF_y])== d[EFF_x]) & (flip_alleles(d[REF_y])== d[REF_x]), d[EFF_x], np.where((flip_alleles(d[EFF_y])== d[REF_x]) & (flip_alleles(d[REF_y]) == d[EFF_x]), d[REF_x], np.nan))))
        d[REF_y]= np.where((d[EFF_y]== d[EFF_x]) & (d[REF_y]== d[REF_x]), d[REF_x], np.where((d[EFF_y]== d[REF_x]) & (d[REF_y]== d[EFF_x]), d[EFF_x], np.where((flip_alleles(d[EFF_y])== d[EFF_x]) & (flip_alleles(d[REF_y])== d[REF_x]), d[REF_x], np.where((flip_alleles(d[EFF_y])== d[REF_x]) & (flip_alleles(d[REF_y]) == d[EFF_x]), d[EFF_x], np.nan))))
        return(d)

d= pd.read_csv(snakemake.input[0], sep= '\t', header= 0)
betas= pd.read_csv(snakemake.input[1], sep= '\t')

d['chr']= d.chr.apply(str)
betas['chr']= betas['chr'].apply(str)
d= pd.merge(betas, d, on= ['chr', 'pos'])
print(d.columns)
d= harmonize_alleles(d, 'ref', 'eff', 'REF', 'EFF')
d.drop_duplicates(['chr', 'pos'], keep=False, inplace= True)
d.drop(['ref', 'eff', 'pos', 'chr', 'EFF', 'REF'], axis= 1, inplace= True)
d= d.iloc[:,1:].multiply(d['beta'], axis="index")

d= pd.DataFrame(d.sum(axis=0))

if 'haplotype' in snakemake.input[0]:
        d['PREG_ID']= d.index
        d.columns.values[0]= snakemake.wildcards.haplo
else:
        d['IID']= d.index
        d.columns.values[0]= 'BW_fetal_effect_GRS'


d.to_csv(snakemake.output[0], sep ='\t', header= True, index= False)

