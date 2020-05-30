#!/usr/bin/python3

import pandas as pd
import numpy as np

coh= snakemake.wildcards.cohort

def pheno_harvest():
	mfr= pd.read_csv(snakemake.input[0], sep= '\t', header= 0)
	trio= pd.read_csv(snakemake.input[1], sep= '\t', header= None, names= ['PREG_ID', 'Child', 'Father', 'Mother'])
	trio.dropna(subset= ['PREG_ID'], inplace= True)
	mfr['PREG_ID']= coh + '_' + mfr.PREG_ID_1724.astype(int).map(str)
#	trio['PREG_ID']= trio.PREG_ID.str.replace('.0', '')
	d= pd.merge(mfr, trio, on= ['PREG_ID'], how= 'inner')
	d= d[(d['FLERFODSEL']==0)]
	d= d[(d['DODKAT']<6) | (d['DODKAT']>10)]
	d= d[(d['SVLEN_UL_DG']< 308)]
	d= d[(d['SVLEN_UL_DG']> 154)]
	d.dropna(subset= ['SVLEN_UL_DG'], inplace= True)
	d= d.sample(frac=1)
	flag= pd.read_csv(snakemake.input[2], sep= '\t', header= 0)
	pca_out= [line.strip() for line in open(snakemake.input[3], 'r')]
	flag= flag[(flag['genotypesOK']== True) & (flag['phenotypesOK']== True)  & (~flag['IID'].isin(pca_out))]
#	d= d.loc[(d.Child.isin(flag.IID)) & (d.Mother.isin(flag.IID)) & (d.Father.isin(flag.IID)), :]
	d['Child']= np.where(d.Child.isin(flag.IID), d.Child, np.nan)
	d['Mother']= np.where(d.Mother.isin(flag.IID), d.Mother, np.nan)
	d['Father']= np.where(d.Father.isin(flag.IID), d.Father, np.nan)
	d['spont']= np.where((d.FSTART==1) | (((d.KSNITT.isnull()) | (d.KSNITT>1)) & ((d.KSNITT_PLANLAGT.isnull()) | (d.KSNITT_PLANLAGT==1)) & (d.INDUKSJON_PROSTAGLANDIN==0) & (d.INDUKSJON_ANNET==0) & (d.INDUKSJON_OXYTOCIN==0) & (d.INDUKSJON_AMNIOTOMI==0)) | (~d.VANNAVGANG.isnull()) , 1, 0)
	d['PROM']= np.where(d.VANNAVGANG.isnull(), 0, 1)
	d['PARITY0']= np.where(d.PARITET_5==0, 1, 0)
	d['cohort']= coh
	d.drop_duplicates(subset= ['Mother'], keep= 'first', inplace= True)
	d['OBint']= np.where(d.spont==0, 1, 0)
	d['MOR_ALDER']= d.MORS_ALDER
	d['FAR_ALDER']= d.FARS_ALDER_KAT_K8
	d= d[['Child', 'Mother', 'Father', 'PREG_ID', 'spont', 'OBint', 'PROM', 'SVLEN_UL_DG', 'PARITY0', 'cohort', 'MOR_ALDER', 'FAR_ALDER']]
	return d

def pheno_rotterdam():
	mfr= pd.read_csv(snakemake.input[0], sep= '\t', header= 0)
	trio= pd.read_csv(snakemake.input[1], sep= '\t', header= None, names= ['PREG_ID', 'Child', 'Father', 'Mother'])
	trio.dropna(subset= ['PREG_ID'], inplace= True)
	mfr['PREG_ID']= coh + '_' + mfr.PREG_ID_315.astype(int).map(str)
	d= pd.merge(mfr, trio, on= ['PREG_ID'], how= 'inner')
	d= d[(d['FLERFODSEL']=='Enkeltfødsel')]
	d= d[d['DODKAT'].str.contains('Levendefødt')]
	d= d[(d['SVLEN_UL_DG']< 308)]
	d= d[(d['SVLEN_UL_DG']> 154)]
	d.dropna(subset= ['SVLEN_UL_DG'], inplace= True)
	flag= pd.read_csv(snakemake.input[2], sep= '\t', header= 0)
	pca_out= [line.strip() for line in open(snakemake.input[3], 'r')]
	flag= flag[(flag['genotypesOK']== True) & (flag['phenoOK']== True) & (~flag['IID'].isin(pca_out))]
	d['Child']= np.where(d.Child.isin(flag.IID), d.Child, np.nan)
	d['Mother']= np.where(d.Mother.isin(flag.IID), d.Mother, np.nan)
	d['Father']= np.where(d.Father.isin(flag.IID), d.Father, np.nan)
	d['spont']= np.where(((d.FSTART=='Spontan') | (d.FSTART== '')) | (((d.KSNITT=='') | (d.KSNITT== 'Uspesifisert') | (d.KSNITT== 'Akutt keisersnitt')) & (d.INDUKSJON_PROSTAGLANDIN=='Nei') & (d.INDUKSJON_ANNET=='Nei') & (d.INDUKSJON_OXYTOCIN=='Nei') & (d.INDUKSJON_AMNIOTOMI=='Nei')) | (~d.VANNAVGANG.isnull()), 1, 0)
	d['PROM']= np.where(d.VANNAVGANG.isnull(), 0, 1)
	d['PARITY0']= np.where(d.PARITET_5=='0 (førstegangsfødende)', 1, 0)
	d['cohort']= coh
	d.drop_duplicates(subset= ['Mother'], keep= 'first', inplace= True)
	d['OBint']= np.where(d.spont==0, 1, 0)
	d['MOR_ALDER']= d.FAAR - d.MOR_FAAR
	d['FAR_ALDER']= pd.Categorical(d.FARS_ALDER_KAT_K8).codes + 1
	d['FAR_ALDER']= np.where(d.FAR_ALDER== 0, np.nan, d.FAR_ALDER)
	d= d[['Child', 'Mother', 'Father', 'PREG_ID', 'spont', 'PROM', 'OBint', 'SVLEN_UL_DG', 'PARITY0', 'cohort', 'MOR_ALDER', 'FAR_ALDER']]
	return d


if 'harvest' in coh:
	d= pheno_harvest()
if 'harvest' not in coh:
	d= pheno_rotterdam()

d.to_csv(snakemake.output[0], sep= '\t', index= False, header= True)


