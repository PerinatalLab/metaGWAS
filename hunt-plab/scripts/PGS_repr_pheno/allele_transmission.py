import pandas as pd
import numpy as np

PREG_ID= 'PREG_ID'
Sentrix= 'IID'

def format_df(df):
        df.columns= d.PREG_ID
        df[['chr', 'pos', 'ref', 'eff']]= varnames.str.split(':', expand= True)
        cols = list(df.columns.values)
        cols= cols[-4:] + cols[:-4]
        df= df[cols]
        return df

def phase_by_transmission(MT_list, PT_list):
        # Credit goes to Alistair Miles - code adapted from scikit-allel

        n_variants = len(varnames)
        n_samples = fets.shape[1]
        n_progeny = 1
        max_allele = 1

        # setup intermediates
        mac = np.zeros(max_allele + 1, dtype='u1')  # maternal allele counts
        pac = np.zeros(max_allele + 1, dtype='u1')  # paternal allele counts

        # setup outputs

        for index, row in d.iterrows():
                mat1= np.array([int(i[0]) for i in moms[row['MOR_PID']]])
                mat2= np.array([int(i[2]) for i in moms[row['MOR_PID']]])
                pat1= np.array([int(i[0]) for i in dads[row['FAR_PID']]])
                pat2= np.array([int(i[2]) for i in dads[row['FAR_PID']]])
                MT= np.array([int(i[0]) for i in fets[row['BARN_PID']]])
                PT= np.array([int(i[2]) for i in fets[row['BARN_PID']]]) 
                for i in range(n_variants):
                        ma1 = mat1[i]  # maternal allele 1
                        ma2 = mat2[i]  # maternal allele 2
                        pa1 = pat1[i]  # paternal allele 1
                        pa2 = pat2[i]  # paternal allele 2
                        # check for any missing calls in parents
                        if ma1 < 0 or ma2 < 0 or pa1 < 0 or pa2 < 0:
                                continue
                        # parental allele counts
                        mac[:] = 0  # reset to zero
                        pac[:] = 0  # reset to zero
                        mac[ma1] = 1
                        mac[ma2] = 1
                        pac[pa1] = 1
                        pac[pa2] = 1
                        a1 = MT[i]
                        a2 = PT[i]
                        if a1 < 0 or a2 < 0:  # child is missing
                                continue
                        elif a1 == a2:  # child is homozygous
                                if mac[a1] > 0 and pac[a1] > 0:  # Mendelian consistent
                                        # trivially phase the child
                                        continue
                        else:  # child is heterozygous
                                if mac[a1] > 0 and pac[a1] == 0 and pac[a2] > 0:
                                # allele 1 is unique to mother, no need to swap
                                        continue
                                elif mac[a2] > 0 and pac[a2] == 0 and pac[a1] > 0:
                                        # allele 2 is unique to mother, swap child alleles
                                        MT[i] = a2
                                        PT[i] = a1
                                elif pac[a1] > 0 and mac[a1] == 0 and mac[a2] > 0:
                                        # allele 1 is unique to father, swap child alleles
                                        MT[i] = a2
                                        PT[i] = a1
                                elif pac[a2] > 0 and mac[a2] == 0 and mac[a1] > 0:
                                        # allele 2 is unique to father, no need to swap
                                        continue
                MT_list.append(pd.Series(MT))
                PT_list.append(pd.Series(PT))
        return MT_list, PT_list

d= pd.read_csv(snakemake.input[2], sep= '\t', header= 0)
d.dropna(axis= 0, inplace= True)

fets= pd.read_csv(snakemake.input[0], sep='\t', header= 0)
moms= pd.read_csv(snakemake.input[1], sep= '\t', header= 0)
dads= pd.read_csv(snakemake.input[3], sep='\t', header= 0)

d['BARN_PID']= d.BARN_PID.apply(int).apply(str)
d['MOR_PID']= d.MOR_PID.apply(int).apply(str)
d['FAR_PID']= d.FAR_PID.apply(int).apply(str)

d= d.loc[d.BARN_PID.isin(fets.columns), :]
d= d.loc[d.MOR_PID.isin(moms.columns), :]
d= d.loc[d.FAR_PID.isin(dads.columns), :]

varnames= fets.chr.map(str) + ':' + fets.pos.map(str) + ':' + fets.ref + ':' + fets.eff

fets= fets[d.BARN_PID]
moms= moms[d.MOR_PID]
dads= dads[d.FAR_PID]


fets.replace('/', '|', inplace= True)
moms.replace('/', '|', inplace= True)
dads.replace('/', '|', inplace= True)

moms1= moms.copy()
dads1= dads.copy()

moms= moms.loc[:,~moms.columns.duplicated()]
dads= dads.loc[:,~dads.columns.duplicated()]

MT_list= list()
PT_list= list()

MT_list, PT_list= phase_by_transmission(MT_list, PT_list)
MT= pd.concat(MT_list, axis= 1)
PT= pd.concat(PT_list, axis= 1)

moms1= np.where(moms1== '0|0', 0, np.where((moms1== '0|1') | (moms1== '1|0'), 1, np.where(moms1== '1|1', 2, np.nan)))
dads1= np.where(dads1== '0|0', 0, np.where((dads1== '0|1') | (dads1== '1|0'), 1, np.where(dads1== '1|1', 2, np.nan)))

MnT= moms1 - MT
PnT= dads1 - PT

d['PREG_ID']= 'PREG' + '_' + d.MOR_PID.apply(str) + '_' + d.BARN_PID.apply(str) + '_' + d.FAR_PID.apply(str) 


MT= format_df(MT)
MnT= format_df(MnT)
PT= format_df(PT)
PnT= format_df(PnT)

MT.to_csv(snakemake.output[0], header= True, sep= '\t', index= False)
MnT.to_csv(snakemake.output[1], header= True, sep= '\t', index= False)
PT.to_csv(snakemake.output[2], header= True, sep= '\t', index= False)
PnT.to_csv(snakemake.output[3], header= True, sep= '\t', index= False)
