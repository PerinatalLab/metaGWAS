import pandas as pd
import numpy as np

d= pd.read_csv(snakemake.input[0], sep= '\t', header=0, compression= 'gzip')
d= d.loc[~d['#chrom'].str.contains('_'), :]
d['a1']= d.alts.str.split(',').str[0]
d['a2']= d.alts.str.split(',').str[1]
d['#chrom']= d['#chrom'].str.replace('chr', '')
d['POS']= np.where(d.ref.str.len() < d.alts.str.len(), d.chromStart, d.chromEnd)
d['ref']= np.where(d.ref.str.len()< d.alts.str.len(), 'I', d.ref)
d['ref']= np.where(d.ref.str.len() > d.alts.str.len(), 'D', d.ref)
d['a1']= np.where(d.ref== 'I', 'D', d.a1)
d['a1']= np.where(d.ref== 'D', 'I', d.a1)
df= d.copy()
df= df.loc[df.a2!= '', :]
d.loc[d.ref > d.a1, ['ref', 'a1']] = d.loc[d.ref > d.a1, ['a1', 'ref']].values

d['ID']= d['#chrom'] + ':' + d['POS'].astype(int).astype(str) + ':' + d.ref + ':' + d.a1
df.loc[df.ref > df.a2, ['ref', 'a2']] = df.loc[df.ref > df.a2, ['a2', 'ref']].values
df['ID']= df['#chrom'] + ':' + df['POS'].astype(int).astype(str) + ':' + df.ref + ':' + df.a2
df= df[['ID', 'name']]
d= d[['ID', 'name']]
d= pd.concat([d, df])

# Read RSIDs from HRC
x= pd.read_csv(snakemake.input[1], sep= '\t', header=0, usecols= ['#CHROM', 'POS', 'ID', 'REF', 'ALT'])
x.columns= ['CHROM', 'POS', 'name', 'REF', 'ALT']
x= x.loc[x.name!= '.', :]

x['CHROM']= np.where(x.CHROM== 'X', '23', x.CHROM)
x['CHROM']= x.CHROM.apply(str)

x.loc[x.REF > x.ALT, ['REF', 'ALT']] = x.loc[x.REF > x.ALT, ['ALT', 'REF']].values
x['ID']= x['CHROM'] + ':' + x['POS'].astype(int).astype(str) + ':' + x.REF + ':' + x.ALT
x= x[['ID', 'name']]
x= x.loc[~x.ID.isin(d.ID), :]

d= pd.concat([d, x])

d.to_csv(snakemake.output[0], sep= '\t', header= True, index= False)

