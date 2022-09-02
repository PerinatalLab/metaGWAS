library(data.table)
library(dplyr)
library(tidyr)
coh= sub('_mfr.csv', '', unlist(strsplit(snakemake@input[[1]], 'pheno/'))[2])
print(coh)
mfr= fread(snakemake@input[[1]])
trio= fread(snakemake@input[[2]], h=F)
names(trio)= c('PREG_ID', 'Child', 'Father', 'Mother')
trio= trio[!is.na(trio$PREG_ID), ]
trio= trio[trio$PREG_ID!= '', ]
flag= fread(snakemake@input[[3]])
out= fread(snakemake@input[[4]])


if (grepl('harvest', snakemake@input[[1]])) {

mfr= filter(mfr, FLERFODSEL==0, (DODKAT<6 | DODKAT>10), is.na(IVF), ABRUPTIOP==0, PLACENTA_PREVIA==0, is.na(PREEKL), EKLAMPSI==0, FOSTERV_POLYHYDRAMNION==0, FOSTERV_OLIGOHYDRAMNION== 0, is.na(DIABETES_MELLITUS), HYPERTENSJON_KRONISK==0, REUM_ARTRITT==0, KSNITT_TIDLIGERE==0, C00_MALF_ALL==0)
mfr= mutate(mfr, spont= as.numeric((FSTART==1 | (is.na(KSNITT) | KSNITT>1) &
                (is.na(KSNITT_PLANLAGT) | KSNITT_PLANLAGT==1)) & (INDUKSJON_PROSTAGLANDIN==0) & (INDUKSJON_ANNET==0) & (INDUKSJON_OXYTOCIN==0) & (is.na(VANNAVGANG))), PROM= as.numeric(!is.na(VANNAVGANG)), RESPIRATORISK_DISTR= as.numeric(RESPIRATORISK_DISTR== 1), INTRAKRANIELL_BLODN= as.numeric(INTRAKRANIELL_BLODN== 1), SYSTEMISKANTIBIOTIKA= as.numeric(SYSTEMISKANTIBIOTIKA== 1), ICTERUS= as.numeric(ICTERUS== 1), OVERFLYTTET= as.numeric(OVERFLYTTET== 2), KJONN= as.numeric(KJONN== 1))

flag= filter(flag, genotypesOK== TRUE, phenotypesOK== TRUE)
flag= filter(flag, !(IID %in% out$V2))
mfr$PREG_ID_1724= paste(coh, mfr$PREG_ID_1724, sep= '_')
mfr= inner_join(mfr, trio, by= c('PREG_ID_1724'= 'PREG_ID'))
mfr$Child= ifelse(mfr$Child %in% flag$IID, mfr$Child, '')
mfr$Mother= ifelse(mfr$Mother %in% flag$IID, mfr$Mother, '')
mfr$Father= ifelse(mfr$Father %in% flag$IID, mfr$Father, '')
#names(mfr)[names(mfr) == 'SentrixID_1']= 'IID'
names(mfr)[names(mfr) == 'PREG_ID_1724']= 'PREG_ID'
mfr$PARITY= as.numeric(mfr$PARITET_5== 0)

}

if (!grepl('harvest', snakemake@input[[1]])) {

mfr= filter(mfr, FLERFODSEL=='Enkeltfødsel', grepl('Levendefødt', DODKAT), is.na(IVF), ABRUPTIOP== 'Nei', PLACENTA_PREVIA=='Nei', is.na(PREEKL), EKLAMPSI=='Nei', FOSTERV_POLYHYDRAMNION=='Nei', FOSTERV_OLIGOHYDRAMNION== 'Nei', is.na(DIABETES_MELLITUS), HYPERTENSJON_KRONISK=='Nei', REUM_ARTRITT=='Nei', C00_MALF_ALL=='Nei')
mfr= mutate(mfr, spont= as.numeric(((FSTART=='Spontan' | FSTART== '') | (KSNITT=='' | KSNITT== 'Uspesifisert' | KSNITT== 'Akutt keisersnitt')) & (INDUKSJON_PROSTAGLANDIN=='Nei') & (INDUKSJON_ANNET=='Nei') & (INDUKSJON_OXYTOCIN=='Nei') & (INDUKSJON_AMNIOTOMI=='Nei')), PROM= as.numeric(!is.na(VANNAVGANG)), RESPIRATORISK_DISTR= as.numeric(RESPIRATORISK_DISTR== 'Ja'), INTRAKRANIELL_BLODN= as.numeric(INTRAKRANIELL_BLODN== 'Ja'), SYSTEMISKANTIBIOTIKA= as.numeric(SYSTEMISKANTIBIOTIKA== 'Ja'), ICTERUS= as.numeric(ICTERUS== 'Ja'), OVERFLYTTET= as.numeric(OVERFLYTTET== 'Ja'), KJONN= as.numeric(KJONN== 'Pike'))

flag= filter(flag, genotypesOK== TRUE, phenoOK== TRUE)
flag= filter(flag, !(IID %in% out$V2))

mfr$PREG_ID_315= paste(coh, mfr$PREG_ID_315, sep= '_')
mfr= inner_join(mfr, trio, by= c('PREG_ID_315'= 'PREG_ID'))
mfr$Child= ifelse(mfr$Child %in% flag$IID, mfr$Child, '')
mfr$Mother= ifelse(mfr$Mother %in% flag$IID, mfr$Mother, '')
mfr$Father= ifelse(mfr$Father %in% flag$IID, mfr$Father, '')

#names(mfr)[names(mfr) == 'SentrixID']= 'IID'
names(mfr)[names(mfr) == 'PREG_ID_315']= 'PREG_ID'
mfr$PARITY= as.numeric(mfr$PARITET_5== '0 (førstegangsfødende)')

}

mfr= filter(mfr, SVLEN_UL_DG>= 140, SVLEN_UL_DG<=310, !is.na(SVLEN_UL_DG))
mfr= mfr[order(mfr$PARITET_5, decreasing= T), ]
#mfr= mfr[!duplicated(mfr$IID), ]
mfr= select(mfr, Mother, Child, Father, PREG_ID, SVLEN_UL_DG, spont, PROM, PARITY, KJONN, MORS_ALDER_KAT_K6, RESPIRATORISK_DISTR, INTRAKRANIELL_BLODN, SYSTEMISKANTIBIOTIKA, ICTERUS, OVERFLYTTET, VEKT, SVLEN_SM_DG)
mfr$cohort= coh

write.table(mfr, snakemake@output[[1]], sep= '\t', quote= FALSE, row.names= FALSE, col.names= TRUE)

