library(data.table)
library(dplyr)
library(tidyr)

mfr= fread(snakemake@input[[1]])
ids= fread(snakemake@input[[2]], h= F)

names(ids)= c('IID', 'BATCH', 'PREG_ID', 'Role')

if (!grepl('TED', snakemake@output[[1]])) {

ids= filter(ids, BATCH != 'TED')

} else {
ids= filter(ids, BATCH == 'TED')
ids= group_by(ids, PREG_ID, Role) %>% filter(row_number()== 1)
}

ids= group_by(ids, PREG_ID, Role) %>% filter(row_number()== 1)

ids= ids[!duplicated(ids[,c('PREG_ID', 'IID')]), ]

flag= fread(snakemake@input[[3]])
flag= filter(flag, genotypesOK== T, phenoOK== T)

#out= readLines(snakemake@input[[4]])


#ids= filter(ids, !(IID %in% out), 
ids= inner_join(ids, flag[, 'IID'], 'IID')
ids= spread(ids, key= Role, value= IID)

mfr= inner_join(mfr, ids, by= c('PREG_ID_1724'= 'PREG_ID'))


mfr= filter(mfr, is.na(FLERFODSEL), is.na(DAAR) | DAAR != FAAR, is.na(ABRUPTIOP), is.na(PLACENTA_PREVIA), is.na(PREEKL), is.na(EKLAMPSI), is.na(DIABETES_MELLITUS), is.na(MISD))

mfr= mutate(mfr, spont= as.numeric(((FSTART==1 | is.na(FSTART)) | (is.na(KSNITT) | KSNITT== 2 | KSNITT== 9)) & (is.na(INDUKSJON_PROSTAGLANDIN)) & (is.na(INDUKSJON_ANNET)) & (is.na(INDUKSJON_OXYTOCIN)) & (is.na(INDUKSJON_AMNIOTOMI))), KJONN= as.numeric(KJONN== 1))

names(mfr)[names(mfr) == 'SentrixID']= 'IID'
names(mfr)[names(mfr) == 'PREG_ID_1724']= 'PREG_ID'
mfr$PARITY= as.numeric(mfr$PARITET_5== 1)

mfr$SVLEN= with(mfr, ifelse(SVLEN_UL_DG< 140 | SVLEN_UL_DG> 310, NA, SVLEN_UL_DG))
mfr$SVLEN= with(mfr, ifelse(is.na(SVLEN_UL_DG), SVLEN_SM_DG, SVLEN_UL_DG))

mfr= filter(mfr, SVLEN>= 140, SVLEN<=310, !is.na(SVLEN))
mfr= arrange(mfr, desc(spont), Mother, Child, Father)

mfr$cohort= mfr$BATCH

trios= fread(snakemake@input[[5]])
mfr= filter(mfr, as.numeric(PREG_ID) %in% as.numeric(trios$PREG_ID))
mfr$MORS_ALDER= mfr$FAAR - mfr$MOR_FAAR 
moms= mfr[!duplicated(mfr$Mother, incomparables= NA), ]
fets= mfr[!duplicated(mfr$Child, incomparables= NA), ]
dads= mfr[!duplicated(mfr$Father, incomparables= NA), ]

moms= select(moms, Mother, PREG_ID, SVLEN, SVLEN_UL_DG, spont,PARITY, KJONN, MORS_ALDER, VEKT, SVLEN_SM_DG, cohort)

fets= select(fets, Child, PREG_ID, SVLEN, SVLEN_UL_DG, spont, PARITY, KJONN, MORS_ALDER, VEKT, SVLEN_SM_DG, cohort)

dads= select(dads, Father, PREG_ID, SVLEN, SVLEN_UL_DG, spont, PARITY, KJONN, MORS_ALDER, VEKT, SVLEN_SM_DG, cohort)

names(moms)[1]= 'IID'
names(fets)[1]= 'IID'
names(dads)[1]= 'IID'

fwrite(moms, snakemake@output[[1]], sep= '\t')
fwrite(fets, snakemake@output[[2]], sep= '\t')
fwrite(dads, snakemake@output[[3]], sep= '\t')

