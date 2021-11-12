import pandas as pd
import numpy as np
import re

#d= pd.read_csv(snakemake.input[0], sep= '\t', header= 0)

#d['Allele1']= d['Allele1'].str.upper()
#d['Allele2']= d['Allele2'].str.upper()
#d= d.loc[(d.TOTALSAMPLESIZE> (d['TOTALSAMPLESIZE'].max())/ 2), :]
#d[['CHR', 'POS', 'REF','EFF', 'SNP']]= d['MarkerName'].str.split(':', expand= True)
#d['CHR']= d['CHR'].astype(str).astype(int)
#d['POS']= d['POS'].astype(str).astype(int)
#d= d[['CHR', 'POS', 'Allele1', 'Allele2', 'TOTALSAMPLESIZE', 'Freq1', 'Effect', 'StdErr', 'P-value']]
#d.columns= ['CHR', 'POS', 'EFF', 'REF', 'TOTALSAMPLESIZE', 'EAF', 'BETA', 'SE', 'pvalue']
#d['BETA']=np.where(d.REF > d.EFF, -1* d.BETA, d.BETA)
#d['EAF']= np.where(d.REF > d.EFF, 1 - d.EAF, d.EAF)

#d['CHR']= d['CHR'].astype(str).astype(int)
#d['POS']= d['POS'].astype(str).astype(int)

#d['pvalue']= d['pvalue'].astype(str).astype(float)

#d.loc[d.REF > d.EFF, ['REF', 'EFF']] = d.loc[d.REF > d.EFF, ['EFF', 'REF']].values
#d['ID']= d.CHR.astype(int).astype(str) + ':' + d.POS.astype(int).astype(str) + ':' + d.REF + ':' + d.EFF

#d= d.loc[((d.pvalue>0) & (d.pvalue <1)), :]

col_list= ['IMPACT', 'DISTANCE', 'SYMBOL', 'SYMBOL_SOURCE', 'BIOTYPE']
df_list= list()

for vep in pd.read_csv(snakemake.input[1], sep= '\t', header= None, names= ['Variation', 'Location', 'Allele', 'Gene', 'Feature', 'Feature_type', 'Consequence', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'Extra'], comment= '#', chunksize= 100000):
	for i in col_list:
		vep[i]= vep['Extra'].apply(lambda y: dict([(x.split('=', 1)) for x in re.split(';(?=\w)', y) if x.find('=') > -1])[i] if i in y else '')
	vep= vep[['Variation', 'Location', 'Existing_variation', 'Gene', 'SYMBOL', 'Consequence', 'IMPACT', 'DISTANCE', 'SYMBOL_SOURCE', 'BIOTYPE']]
	vep.columns= ['ID', 'Location', 'RSID', 'Gene', 'SYMBOL', 'Consequence', 'IMPACT', 'DISTANCE', 'SYMBOL_SOURCE', 'BIOTYPE']
	vep['BIOTYPE1']= np.where(vep.BIOTYPE== 'protein_coding', 0, np.where(vep.BIOTYPE.str.contains('pseudo'), 2, 1))
	vep['DISTANCE']= np.where(vep.DISTANCE== '', 0, vep.DISTANCE)
	vep[['chr', 'pos', 'All']]= vep.ID.str.split('_', expand= True)
	vep[['EFF', 'REF']]= vep.All.str.split('/', expand= True)
	vep.loc[vep.REF > vep.EFF, ['REF', 'EFF']] = vep.loc[vep.REF > vep.EFF, ['EFF', 'REF']].values
	vep[['CHR', 'POS']]= vep['Location'].str.split(':', expand= True)
	vep['CHR']= np.where(vep['CHR']== 'X', '23', vep['CHR'])
	vep['ID']= vep.CHR.astype(int).astype(str) + ':' + vep.POS.astype(int).astype(str) + ':' + vep.REF + ':' + vep.EFF
	vep= vep[['ID', 'RSID', 'Gene', 'SYMBOL', 'Consequence', 'IMPACT', 'DISTANCE', 'BIOTYPE', 'BIOTYPE1']]
	vep.sort_values(by= ['BIOTYPE1'], ascending= True, inplace= True)
	vep.drop_duplicates(subset= ['ID'], keep= 'first', inplace= True)
	df_list.append(vep)

vep= pd.concat(df_list)

vep.sort_values(by= ['BIOTYPE1'], ascending= True, inplace= True)
vep.drop_duplicates(subset= ['ID'], keep= 'first', inplace= True)
vep= vep[['ID', 'RSID', 'Gene', 'SYMBOL', 'Consequence', 'IMPACT', 'DISTANCE', 'BIOTYPE']]


d= pd.read_csv(snakemake.input[0], sep= '\t', header= 0)
d['Allele1']= d['Allele1'].str.upper()
d['Allele2']= d['Allele2'].str.upper()
d= d.loc[(d.TOTALSAMPLESIZE> (d['TOTALSAMPLESIZE'].max())/ 2), :]
d= d.loc[d.TOTALSAMPLESIZE> 66106, :]
d[['CHR', 'POS', 'REF','EFF']]= d['MarkerName'].str.split(':', expand= True)
d['CHR']= d['CHR'].astype(str).astype(int)
d['POS']= d['POS'].astype(str).astype(int)
d= d[['CHR', 'POS', 'Allele1', 'Allele2', 'TOTALSAMPLESIZE', 'Freq1', 'P-value']]
d.columns= ['CHR', 'POS', 'EFF', 'REF', 'TOTALSAMPLESIZE', 'EAF', 'pvalue']
d['EAF']= np.where(d.REF > d.EFF, 1 - d.EAF, d.EAF)
d['CHR']= d['CHR'].astype(str).astype(int)
d['POS']= d['POS'].astype(str).astype(int)
d['pvalue']= d['pvalue'].astype(str).astype(float)
d.loc[d.REF > d.EFF, ['REF', 'EFF']] = d.loc[d.REF > d.EFF, ['EFF', 'REF']].values
d['ID']= d.CHR.astype(int).astype(str) + ':' + d.POS.astype(int).astype(str) + ':' + d.REF + ':' + d.EFF
d= d.loc[((d.pvalue>0) & (d.pvalue <1)), :]
d['MAF']= np.where(d.EAF>0.5, 1 - d.EAF, d.EAF)
d= d.loc[d.MAF>= 0.1, :]
d= pd.merge(d, vep, on= ['ID'], how= 'left')
d.to_csv(snakemake.output[0], header=True, index= False, sep= '\t')
