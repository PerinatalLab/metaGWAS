library(data.table)
library(dplyr)
library(tidyr)

mfr= fread(snakemake@input[[1]])
link= fread(snakemake@input[[2]])
flag= fread(snakemake@input[[3]])
out= fread(snakemake@input[[4]])

link= filter(link, Role== 'Mother')

if (grepl('harvest', snakemake@input[[1]])) {

mfr= filter(mfr, FLERFODSEL==0, (DODKAT<6 | DODKAT>10), is.na(IVF), ABRUPTIOP==0, PLACENTA_PREVIA==0, is.na(PREEKL), EKLAMPSI==0, FOSTERV_POLYHYDRAMNION==0, FOSTERV_OLIGOHYDRAMNION== 0, is.na(DIABETES_MELLITUS), HYPERTENSJON_KRONISK==0, REUM_ARTRITT==0, KSNITT_TIDLIGERE==0, C00_MALF_ALL==0)
mfr= mutate(mfr, spont= as.numeric(FSTART==1 & (is.na(KSNITT) | KSNITT>1) &
                (is.na(KSNITT_PLANLAGT) | KSNITT_PLANLAGT==1)))
mfr= inner_join(mfr, link, by= 'PREG_ID_1724')
flag= filter(flag, genotypesOK== TRUE, phenotypesOK== TRUE)
mfr= filter(mfr, SentrixID_1 %in% flag$IID)
mfr= filter(mfr, !(SentrixID_1 %in% out$V2))
names(mfr)[names(mfr) == 'SentrixID_1']= 'IID'
names(mfr)[names(mfr) == 'PREG_ID_1724']= 'PREG_ID'

}

if (!grepl('harvest', snakemake@input[[1]])) {

mfr= filter(mfr, FLERFODSEL=='Enkeltfødsel', grepl('Levendefødt', DODKAT), is.na(IVF), ABRUPTIOP== 'Nei', PLACENTA_PREVIA=='Nei', is.na(PREEKL), EKLAMPSI=='Nei', FOSTERV_POLYHYDRAMNION=='Nei', FOSTERV_OLIGOHYDRAMNION== 'Nei', is.na(DIABETES_MELLITUS), HYPERTENSJON_KRONISK=='Nei', REUM_ARTRITT=='Nei', C00_MALF_ALL=='Nei')
mfr= mutate(mfr, spont= as.numeric(FSTART=='Spontan' | FSTART== '' & (KSNITT=='' | KSNITT== 'Uspesifisert' | KSNITT== 'Akutt keisersnitt')))

mfr= inner_join(mfr, link, by= 'PREG_ID_315')
flag= filter(flag, genotypesOK== TRUE, phenoOK== TRUE)
mfr= filter(mfr, SentrixID %in% flag$IID)
mfr= filter(mfr, !(SentrixID %in% out$V2))
names(mfr)[names(mfr) == 'SentrixID']= 'IID'
names(mfr)[names(mfr) == 'PREG_ID_315']= 'PREG_ID'

}

p1= mfr
p1$allPTD= ifelse((p1$SVLEN_UL_DG<259 & p1$spont==1), 1, NA)
p1$allPTD= ifelse((p1$SVLEN_UL_DG>= 273 & p1$SVLEN_UL_DG< 294), 0, p1$allPTD)

p1= filter(p1, !is.na(allPTD))
p1= p1[order(p1$PARITET_5, decreasing= F), ]
p1= p1[!duplicated(p1$IID), ]
p1= select(p1, IID, allPTD)

p2= mfr
p2$earlyPTD= ifelse((p2$SVLEN_UL_DG<238 & p2$spont==1), 1, NA)
p2$earlyPTD= ifelse((p2$SVLEN_UL_DG>= 273 & p2$SVLEN_UL_DG< 294), 0, p2$earlyPTD)
p2= filter(p2, !is.na(earlyPTD))
p2= p2[order(p2$PARITET_5, decreasing= F), ]
p2= p2[!duplicated(p2$IID), ]
p2= select(p2, IID, earlyPTD)

p3= mfr
p3$postTerm= ifelse((p3$SVLEN_UL_DG>=273 & p3$SVLEN_UL_DG < 294  & p3$spont==1), 0, NA)
p3$postTerm= ifelse(p3$SVLEN_UL_DG>= 294, 1, p3$postTerm)
p3= filter(p3, !is.na(postTerm))
p3= p3[order(p3$PARITET_5, decreasing= F), ]
p3= p3[!duplicated(p3$IID), ]
p3= select(p3, IID, postTerm)


p4= mfr
p4= filter(p4, spont== 1, SVLEN_UL_DG>= 140, SVLEN_UL_DG<=310)
p4$GAraw= p4$SVLEN_UL_DG
p4= filter(p4, !is.na(GAraw))
p4= p4[order(p4$PARITET_5, decreasing= F), ]
p4= p4[!duplicated(p4$IID), ]
p4$GAnrm= qnorm((rank(p4$GAraw, na.last="keep") - 0.375) / (sum(!is.na(p4$GAraw)) + 0.25))
p4= select(p4, IID, GAraw, GAnrm)

d= full_join(p1, p2, by= 'IID') %>% full_join(., p3, by= 'IID') %>% full_join(., p4, by= 'IID')

d$cohort= snakemake@wildcards[['cohort']]

d[is.na(d)]= 'nan'

ids= fread(snakemake@input[[5]], h=F)
names(ids)= c('FID', 'IID')

d= inner_join(d, ids, by= 'IID')

write.table(d, snakemake@output[[1]], sep= '\t', quote= FALSE, row.names= FALSE, col.names= TRUE)

write.table(d[, c('FID', 'IID')], snakemake@output[[2]], sep= '\t', quote= F, row.names=F, col.names=F)
#write.table(earlyPTD, snakemake@output[[2]], sep= '\t', quote= FALSE, row.names= FALSE, col.names= TRUE)
#write.table(postTerm, snakemake@output[[3]], sep= '\t', quote= FALSE, row.names= FALSE, col.names= TRUE)
#write.table(GAraw, snakemake@output[[4]], sep= '\t', quote= FALSE, row.names= FALSE, col.names= TRUE)
#write.table(GAnrm, snakemake@output[[5]], sep= '\t', quote= FALSE, row.names= FALSE, col.names= TRUE)

