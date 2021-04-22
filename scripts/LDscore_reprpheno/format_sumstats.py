import pandas as pd
import numpy as np
from scipy.special import chdtri
import gzip
import csv

def not_number(s):
	if s != None:
		try:
			float(s)
			return False
		except ValueError:
			return True
	else:
		return True


def select_format(repr_pheno, row):
	'For each wildcard assign the correct formating function.'
	if repr_pheno== 'Preeclampsia':
		rsid, CHR, POS, EAF, N, REF, EFF, BETA, SE, pvalue= preeclampsia(row)
	if repr_pheno== 'POP': 
		rsid, CHR, POS, EAF, N, REF, EFF, BETA, SE, pvalue= POP(row)
	if repr_pheno== 'miscarriage':
		rsid, CHR, POS, EAF, N, REF, EFF, BETA, SE, pvalue= miscarriage(row)
	if repr_pheno== 'GA_fetal':
		rsid, CHR, POS, EAF, N, REF, EFF, BETA, SE, pvalue= fet_GA(row)
	if repr_pheno== 'BW_maternal':
		rsid, CHR, POS, EAF, N, REF, EFF, BETA, SE, pvalue= BW_maternal(row)
	if repr_pheno== 'BW_fetal':
		rsid, CHR, POS, EAF, N, REF, EFF, BETA, SE, pvalue= BW_fetal(row)
	if repr_pheno== 'BW_maternal_effect':
		rsid, CHR, POS, EAF, N, REF, EFF, BETA, SE, pvalue= BW_maternal_adjusted_effect(row)
	if repr_pheno== 'BW_fetal_effect':
		rsid, CHR, POS, EAF, N, REF, EFF, BETA, SE, pvalue= BW_fetal_adjusted_effect(row)
	if repr_pheno== 'leiomyoma_uterus':
		rsid, CHR, POS, EAF, N, REF, EFF, BETA, SE, pvalue= leiomyoma_uterus(row)
	if repr_pheno in ['Oestradiol_fem', 'NLB', 'AFB', 'AMenarche', 'AMenopause', 'endometriosis']:
		rsid, CHR, POS, EAF, N, REF, EFF, BETA, SE, pvalue= UKBB_traits(row)
	if repr_pheno in ['SHBG_fem', 'Testosterone_fem', 'Testosterone_male', 'SHBG_male', 'CBAT_fem', 'CBAT_male']:
		rsid, CHR, POS, EAF, N, REF, EFF, BETA, SE, pvalue= pritchard(row)
	if repr_pheno == 'PCOS':
		rsid, CHR, POS, EAF, N, REF, EFF, BETA, SE, pvalue= PCOS(row)
	return [rsid, CHR, POS, EAF, N, REF, EFF, BETA, SE, pvalue]


def pritchard(row):
	''
	EAF= float(row['A1_FREQ'])
	CHR= row['#CHROM']
	if CHR== 'X': CHR= 23
	if not_number(CHR): return [0, 0, 0 , 0, 0, 0, 0, 0, 0, 0]
	POS= int(row['POS'])
	CHR= int(CHR)
	REF= row['REF']
	EFF= row['ALT']
	N= int(row['OBS_CT'])
	if not_number(row['BETA']): return [0, 0, 0 , 0, 0, 0, 0, 0, 0, 0]
	if not_number(row['SE']): return [0, 0, 0 , 0, 0, 0, 0, 0, 0, 0]
	if not_number(row['P']): return [0, 0, 0 , 0, 0, 0, 0, 0, 0, 0]
	BETA= float(row['BETA'])
	SE= float(row['SE'])
	pvalue= float(row['P'])
	rsid= row['ID']
	return [rsid, CHR, POS, EAF, N, REF, EFF, BETA, SE, pvalue]

def leiomyoma_uterus(row):
	''
	EAF= float(row['EAF'])
	CHR= row['CHR']
	if CHR== 'X': CHR= 23
	if not_number(CHR): return [0, 0, 0 , 0, 0, 0, 0, 0, 0, 0]
	POS= int(row['POS'])
	CHR= int(CHR)
	REF= row['REF']
	EFF= row['EFF']
	N= row['TOTALSAMPLESIZE']
	BETA= float(row['beta'])
	SE= float(row['se'])
	pvalue= float(row['pvalue'])
	rsid= ''
	return [rsid, CHR, POS, EAF, N, REF, EFF, BETA, SE, pvalue]

def preeclampsia(row):
	''
	CHR= row['CHR']
	if CHR== 'X': CHR= 23
	if not_number(CHR): return [0, 0, 0 , 0, 0, 0, 0, 0, 0, 0]
	POS= int(row['POS'])
	CHR= int(CHR)
	REF= row['REF'].upper()
	EFF= row['EFF'].upper()
	N= 4630 + 373345
	rsid= row['rsid']
	BETA= float(row['beta'])
	SE= float(row['se'])
	EAF= float(row['EAF'])
	pvalue= float(row['pvalue'])
	return [rsid, CHR, POS, EAF, N, REF, EFF, BETA, SE, pvalue]

def BW_fetal_adjusted_effect(row):
	'Define each header for Birth weight fetal effect.'
	EAF= float(row['eaf'])
	CHR= row['chr']
	if CHR== 'X': CHR= 23
	CHR= int(CHR)
	POS= int(row['pos'])
	REF= row['nea'].upper()
	if REF== 'R': REF= 'D'
	EFF= row['ea'].upper()
	if EFF== 'R': EFF= 'D'
	BETA= float(row['beta'])
	pvalue= float(row['p'])
	SE= float(row['se'])
	N= int(row['n_ownBW'])
	rsid= row['RSID']
	return [rsid, CHR, POS, EAF, N, REF, EFF, BETA, SE, pvalue]

def BW_maternal_adjusted_effect(row):
	'Define each header for Birth weight fetal effect.'
	EAF= float(row['eaf'])
	CHR= row['chr']
	if CHR== 'X': CHR= 23
	CHR= int(CHR)
	POS= int(row['pos'])
	REF= row['nea'].upper()
	if REF== 'R': REF= 'D'
	EFF= row['ea'].upper()
	if EFF== 'R': EFF= 'D'
	BETA= float(row['beta'])
	pvalue= float(row['p'])
	SE= float(row['se'])
	N= int(row['n_offBW'])
	rsid= row['RSID']
	return [rsid, CHR, POS, EAF, N, REF, EFF, BETA, SE, pvalue]


def BW_maternal(row):
	'Define each header for Birth weight maternal effect.'
	EAF= float(row['eaf'])
	CHR= row['chr']
	if CHR== 'X': CHR= 23
	CHR= int(CHR)
	POS= int(row['pos'])
	REF= row['nea']
	EFF= row['ea']
	if REF== 'R': REF= 'D'
	if EFF== 'R': EFF= 'D'
	BETA= float(row['beta'])
	pvalue= float(row['p'])
	SE= float(row['se'])
	N= int(row['n'])
	rsid= row['SNP']
	return [rsid, CHR, POS, EAF, N, REF, EFF, BETA, SE, pvalue]

def BW_fetal(row):
	'Define each header for Birth weight maternal effect.'
	EAF= float(row['eaf'])
	CHR= row['chr']
	if CHR== 'X': CHR= 23
	CHR= int(CHR)
	POS= int(row['pos'])
	REF= row['nea']
	EFF= row['ea']
	if REF== 'R': REF= 'D'
	if EFF== 'R': EFF= 'D'
	BETA= float(row['beta'])
	pvalue= float(row['p'])
	SE= float(row['se'])
	N= int(row['n'])
	rsid= row['rsid']
	return [rsid, CHR, POS, EAF, N, REF, EFF, BETA, SE, pvalue]


def PCOS(row):
	'Define each header for PCOS excluding 23andme.'
	EAF= float(row['EAF'])
	CHR= row['CHR']
	if CHR== 'X': CHR= 23
	CHR= int(CHR)
	POS= int(row['POS'])
	REF= row['REF']
	EFF= row['EFF']
	BETA= float(row['beta'])
	pvalue= float(row['pvalue'])
	SE= float(row['se'])
	N= int(round(float(row['TOTALSAMPLESIZE'])))
	rsid= ''
	return [rsid, CHR, POS, EAF, N, REF, EFF, BETA, SE, pvalue]

def UKBB_traits(row):
	'Define each header for UKBB traits (hormones).'
	if row['low_confidence_variant']== 'true': return [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
	CHR= row['variant'].split(':')[0]
	if CHR== 'X': CHR= 23
	POS= row['variant'].split(':')[1]
	if any([not_number(t) for t in [row['minor_AF'], CHR, POS, row['beta'], row['pval'], row['se'], row['n_complete_samples']]]): return [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
	CHR= int(CHR)
	POS= int(POS)
	REF= row['variant'].split(':')[2]
	EFF= row['variant'].split(':')[3]
	BETA= float(row['beta'])
	pvalue= float(row['pval'])
	SE= float(row['se'])
	N= int(row['n_complete_samples'])
	if row['minor_allele']== EFF:
		EAF= float(row['minor_AF'])
	else:
		EAF= 1- float(row['minor_AF'])
	rsid= ''
	return [rsid, CHR, POS, EAF, N, REF, EFF, BETA, SE, pvalue]

def AP_repr(row):
	'Define each header for BOLT-LMM sumstats.'
	EAF= float(row['EAF'])
	CHR= row['CHR']
	if CHR== 'X': CHR= 23
	CHR= int(CHR)
	POS= int(row['POS'])
	REF= row['A2']
	EFF= row['A1']
	BETA= float(row['Beta'])
	pvalue= float(row['P'])
	SE= float(row['se'])
	N= row['N']
	rsid= row['SNP']
	return [rsid, CHR, POS, EAF, N, REF, EFF, BETA, SE, pvalue]


def POP(row):
	'Define each header for pelvic organ prolapse.'
	if not row['CHR'].isdigit(): return [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
	EAF= float(row['EAF'])
	MAF= np.where(EAF> 0.5, 1 - EAF, EAF)
	if MAF < 0.005: return [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
	if row['CHR']== 'X': row['CHR']= 23
	CHR= int(row['CHR'])
	POS= int(row['POS'])
	REF= row['REF']
	EFF= row['EFF']
	BETA= float(row['BETA'])
	pvalue= float(row['pvalue'])
	SE= float(row['SE'])
	N= float(row['N'])
	rsid= ''
	return [rsid, CHR, POS, EAF, N, REF, EFF, BETA, SE, pvalue]

def fet_GA(row):
	'Define each header for Fetal gestational duration.'
	EAF= ''
	if row['Chr']== 'X': row['Chr']= 23
	CHR= int(row['Chr'])
	POS= int(row['Pos'])
	REF= row['Non_effect_allele'].upper()
	EFF= row['Effect_allele'].upper()
	BETA= float(row['Effect'])
	pvalue= float(row['P'])
	SE= float(row['StdErr'])
	N= int(row['N'])
	rsid= row['Rsid']
	return [rsid, CHR, POS, EAF, N, REF, EFF, BETA, SE, pvalue]

def miscarriage(row):
	'Define each header for Miscarriage.'
	EAF= row['Freq1']
	CHR= row['MarkerName'].split(':')[0]
	if CHR== 'X': CHR= 23
	CHR= int(CHR)
	POS= int(row['MarkerName'].split(':')[1])
	REF= row['Allele2'].upper()
	EFF= row['Allele1'].upper()
	BETA= float(row['Effect'])
	pvalue= float(row['P-value'])
	SE= float(row['StdErr'])
	N= 49996 + 174109
	rsid= ''
	return [rsid, CHR, POS, EAF, N, REF, EFF, BETA, SE, pvalue]


def format_list(input, output):
	with gzip.open(input, 'rt', newline='') as f:
		print(input)
		dialect = csv.Sniffer().sniff(f.readline(), delimiters= ' \t')
		f.seek(0)
		input_file= csv.DictReader(f, dialect= dialect)
		df_list= list()
		with open(output, 'w') as csvfile:
			writer = csv.writer(csvfile, delimiter= '\t')
			writer.writerow([g for g in ['ID', 'rsid', 'CHR', 'POS', 'EAF', 'N', 'REF', 'EFF', 'BETA', 'SE', 'pvalue']])
		for row in input_file:
			rsid, CHR, POS, EAF, N, REF, EFF, BETA, SE, pvalue= select_format(snakemake.wildcards.repr_pheno, row)
			if CHR== 0: continue
			if len(REF) >1: REF= 'I'
			if len(EFF) >1: EFF= 'I'
			if REF== 'I': EFF= 'D'
			if EFF== 'I': REF= 'D'
			if REF> EFF:
				ID= str(CHR) + ':' + str(POS) + ':' + EFF + ':' + REF
				BETA= -1 * float(BETA)
				ref= EFF
				eff= REF
				EAF= 1 - float(EAF)
			else:
				ID= str(CHR) + ':' + str(POS) + ':' + REF + ':' + EFF
				BETA= float(BETA)
				eff= EFF
				ref= REF
			df_list.append([ID, rsid, CHR, POS, EAF, N, ref, eff, BETA, SE, pvalue])
			if len(df_list)== 1000:
				with open(output, 'a', newline= '') as file_handler:
					writer1= csv.writer(file_handler, delimiter= '\t')
					for item in df_list:
						writer1.writerow(item)
				df_list= list()
	with open(output, 'a', newline= '') as file_handler:
			writer1= csv.writer(file_handler, delimiter= '\t')
			for item in df_list:
				writer1.writerow(item)


format_list(snakemake.input[0], snakemake.output[0])
