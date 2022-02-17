import pandas as pd
import numpy as np

PREG_ID= 'PREG_ID'
Sentrix= 'IID'

def format_df(df):
	df[['chr', 'pos', 'ref', 'eff']]= df['index'].str.split(':', expand= True)
	cols = list(df.columns.values)
	cols= cols[-4:] + cols[:-4]
	df= df[cols]
	df.drop(['index'], axis= 1, inplace= True)
	return df


d= pd.read_csv(snakemake.input[2], sep= '\t', header= 0)
d.dropna(axis= 0, inplace= True)

fets= pd.read_csv(snakemake.input[0], sep='\t', header= 0)
moms= pd.read_csv(snakemake.input[1], sep= '\t', header= 0)
dads= pd.read_csv(snakemake.input[3], sep='\t', header= 0)
d= d.loc[d.Child.isin(fets.columns), :]
d= d.loc[d.Mother.isin(moms.columns), :]
d= d.loc[d.Father.isin(dads.columns), :]

varnames= fets.chr.map(str) + ':' + fets.pos.map(str) + ':' + fets.ref + ':' + fets.eff
fets= fets[d.Child]

h3= np.where((fets== '0|0') | (fets== '0|1'), 0, np.where((fets== '1|0') | (fets== '1|1'), 1, np.nan))
h1= np.where((fets== '0|0') | (fets== '1|0'), 0, np.where((fets== '0|1') | (fets== '1|1'), 1, np.nan))

varnames_moms= moms.chr.map(str) + ':' + moms.pos.map(str) + ':' + moms.ref + ':' + moms.eff

moms= moms[d.Mother]
moms= np.where(moms== '0|0', 0, np.where((moms== '0|1') | (moms== '1|0'), 1, np.where(moms== '1|1', 2, np.nan)))

h2= moms - h1

dads= dads[d.Father]
dads= np.where(dads== '0|0', 0, np.where((dads== '0|1') | (dads== '1|0'), 1, np.where(dads== '1|1', 2, np.nan)))

h4= dads - h3

h1= pd.DataFrame(data= h1, columns= d[PREG_ID], index= varnames).reset_index()
h2= pd.DataFrame(data= h2, columns= d[PREG_ID], index= varnames).reset_index()
h3= pd.DataFrame(data= h3, columns= d[PREG_ID], index= varnames).reset_index()
h4= pd.DataFrame(data= h4, columns= d[PREG_ID], index= varnames).reset_index()

h1= format_df(h1)
h2= format_df(h2)
h3= format_df(h3)
h4= format_df(h4)

h1.to_csv(snakemake.output[0], header= True, sep= '\t', index= False)
h2.to_csv(snakemake.output[1], header= True, sep= '\t', index= False)
h3.to_csv(snakemake.output[2], header= True, sep= '\t', index= False)
h4.to_csv(snakemake.output[3], header= True, sep= '\t', index= False)

