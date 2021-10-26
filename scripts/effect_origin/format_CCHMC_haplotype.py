import pandas as pd
import numpy as np


def flip_beta(df):
        'Flip EFF and REF allele if REF> EFF. Flip beta direction with same condition. Assumed column names: beta, REF, EFF.'
        df['BETA']= np.where(df.REF>df.EFF, -1 * df.BETA, df.BETA)
        df['REF'], df['EFF']= np.where(df.REF> df.EFF, [df.EFF, df.REF], [df.REF, df.EFF])
        return df

def add_ID(x):
	x['REF']= np.where(x.REF.str.len() > x.EFF.str.len(), 'I', x.REF)
	x['EFF']= np.where(x.REF.str.len() < x.EFF.str.len(), 'I', x.EFF)
	x['REF']= np.where(x.EFF== 'I', 'D', x.REF)
	x['EFF']= np.where(x.REF== 'I', 'D', x.EFF)
	x['ID']= np.where(x.REF> x.EFF, x.CHR.apply(str) + ':' + x.POS.apply(str) + ':' + x.EFF + ':' + x.REF, x.CHR.apply(str) + ':' + x.POS.apply(str) + ':' + x.REF + ':' + x.EFF)
	x= flip_beta(x)
	return x


def format_df(x, reg):
	d= pd.read_csv(x, sep= ',', header= 0)
	d['chr']= d.chr.apply(str)
	d= pd.merge(d, reg, left_on= 'chr', right_on= 'CHR')
	d= d.loc[((d.pos >= d.pos1) & (d.pos<= d.pos2)), :]
	h1= d.loc[:, ['chr', 'pos', 'ref', 'alt', 'h1.coef', 'h1.se', 'h1.pval']]
	h1.columns= ['CHR', 'POS', 'REF', 'EFF', 'BETA', 'SE', 'pvalue']
	h2= d.loc[:, ['chr', 'pos', 'ref', 'alt', 'h2.coef', 'h2.se', 'h2.pval']]
	h2.columns= ['CHR', 'POS', 'REF', 'EFF', 'BETA', 'SE', 'pvalue']
	h3= d.loc[:, ['chr', 'pos', 'ref', 'alt', 'h3.coef', 'h3.se', 'h3.pval']]
	h3.columns= ['CHR', 'POS', 'REF', 'EFF', 'BETA', 'SE', 'pvalue']
	h1= add_ID(h1)
	h2= add_ID(h2)
	h3= add_ID(h3)
	h1.to_csv(snakemake.output[0], sep= '\t', header= True, index= False, columns= ['ID', 'CHR', 'POS', 'REF', 'EFF', 'BETA', 'SE', 'pvalue'])
	h2.to_csv(snakemake.output[1], sep= '\t', header= True, index= False, columns= ['ID', 'CHR', 'POS', 'REF', 'EFF', 'BETA', 'SE', 'pvalue'])
	h3.to_csv(snakemake.output[2], sep= '\t', header= True, index= False, columns= ['ID', 'CHR', 'POS', 'REF', 'EFF', 'BETA', 'SE', 'pvalue'])
	print('Completed file:' + x)

regions= pd.read_csv(snakemake.input[0], sep= '\t', header= 0)

format_df(snakemake.input[1], regions)
