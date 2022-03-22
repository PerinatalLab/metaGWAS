library(data.table)
library(dplyr)
library(tidyr)

mfr= fread(snakemake@input[[1]])
ids= fread(snakemake@input[[2]], h= F)

names(ids)= c('PREG_ID', 'IID', 'BATCH', 'Role')

if (!grepl('TED', snakemake@output[[1]])) {

ids= filter(ids, BATCH != 'TED')

} else {
ids= filter(ids, BATCH == 'TED')
ids= group_by(ids, PREG_ID, Role) %>% filter(row_number()== 1)
}

ids= ids[!duplicated(ids[,c('PREG_ID', 'IID')]), ]

flag= fread(snakemake@input[[3]])
flag= filter(flag, genotypesOK== T, phenoOK== T)

out= readLines(snakemake@input[[4]])


ids= filter(ids, !(IID %in% out), IID %in% flag$IID)

ids= spread(ids, key= Role, value= IID)

mfr= inner_join(mfr, ids, by= c('PREG_ID_315'= 'PREG_ID'))

mfr= filter(mfr, FLERFODSEL=='Enkeltfødsel', grepl('Levendefødt', DODKAT), is.na(IVF), ABRUPTIOP== 'Nei', PLACENTA_PREVIA=='Nei', is.na(PREEKL), EKLAMPSI=='Nei', FOSTERV_POLYHYDRAMNION=='Nei', FOSTERV_OLIGOHYDRAMNION== 'Nei', is.na(DIABETES_MELLITUS), HYPERTENSJON_KRONISK=='Nei', REUM_ARTRITT=='Nei', C00_MALF_ALL=='Nei')
mfr= mutate(mfr, spont= as.numeric(((FSTART=='Spontan' | FSTART== '') | (KSNITT=='' | KSNITT== 'Uspesifisert' | KSNITT== 'Akutt keisersnitt')) & (INDUKSJON_PROSTAGLANDIN=='Nei') & (INDUKSJON_ANNET=='Nei') & (INDUKSJON_OXYTOCIN=='Nei') & (INDUKSJON_AMNIOTOMI=='Nei')), PROM= as.numeric(!is.na(VANNAVGANG)), RESPIRATORISK_DISTR= as.numeric(RESPIRATORISK_DISTR== 'Ja'), INTRAKRANIELL_BLODN= as.numeric(INTRAKRANIELL_BLODN== 'Ja'), SYSTEMISKANTIBIOTIKA= as.numeric(SYSTEMISKANTIBIOTIKA== 'Ja'), ICTERUS= as.numeric(ICTERUS== 'Ja'), OVERFLYTTET= as.numeric(OVERFLYTTET== 'Ja'), KJONN= as.numeric(KJONN== 'Pike'))

names(mfr)[names(mfr) == 'SentrixID']= 'IID'
names(mfr)[names(mfr) == 'PREG_ID_315']= 'PREG_ID'
mfr$PARITY= as.numeric(mfr$PARITET_5== '0 (førstegangsfødende)')

mfr$SVLEN= with(mfr, ifelse(is.na(SVLEN_UL_DG), SVLEN_SM_DG, SVLEN_UL_DG))

mfr= filter(mfr, SVLEN_UL_DG>= 140, SVLEN_UL_DG<=310, !is.na(SVLEN))
mfr= arrange(mfr, desc(spont), Mother, Child, Father)

mfr= mfr[!duplicated(mfr$Mother, incomparables= NA), ]
mfr= mfr[!duplicated(mfr$Child, incomparables= NA), ]
mfr= mfr[!duplicated(mfr$Father, incomparables= NA), ]

mfr$cohort= mfr$BATCH

mfr= select(mfr, Mother, Child, Father, PREG_ID, SVLEN, SVLEN_UL_DG, spont, PROM, PARITY, KJONN, MORS_ALDER_KAT_K6, RESPIRATORISK_DISTR, INTRAKRANIELL_BLODN, SYSTEMISKANTIBIOTIKA, ICTERUS, OVERFLYTTET, VEKT, SVLEN_SM_DG, cohort)

fwrite(mfr, snakemake@output[[1]], sep= '\t')
