import pandas as pd
import numpy as np

d= pd.read_csv(snakemake.input[0], header= 0, delim_whitespace= True)

gest_age_pos= 123112292
bw_pos= 123065778

df= d.copy()

df.columns= ['CHR_B', 'BP_B', 'SNP_B', 'MAF_B', 'CHR_A', 'BP_A', 'SNP_A', 'MAF_A', 'R2']

d= pd.concat([d, df])

d= d.loc[d.R2>= 0.8, :]

pos_list= [pos for pos in d.loc[d.BP_A== gest_age_pos, :]['BP_B']] + [gest_age_pos]

pos_list2= [pos for pos in d.loc[d.BP_A== bw_pos, :]['BP_B']] + [bw_pos]


df_list= list()
for i in pos_list:
	df_list.append([pos for pos in d.loc[d.BP_A== i]['BP_B']] + [i])

gest_age_list= set.intersection(*map(set,df_list))


df_list= list()
for i in pos_list2:
        df_list.append([pos for pos in d.loc[d.BP_A== i]['BP_B']] + [i])

bw_list= set.intersection(*map(set,df_list))

d['haplotype']= np.where((d.BP_A.isin(gest_age_list)) & (d.BP_B.isin(gest_age_list)), 'Gestational duration', None)
d['haplotype']= np.where((d.BP_A.isin(bw_list)) & (d.BP_B.isin(bw_list)), 'Birth weight', d.haplotype)

d.dropna(subset= ['haplotype'], inplace= True)

d.to_csv(snakemake.output[0], sep= '\t', header= True, index= False)
