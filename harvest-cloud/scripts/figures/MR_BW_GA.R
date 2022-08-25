library(dplyr)
library(data.table)
library(MendelianRandomization)

d= fread(snakemake@input[[1]])
d= filter(d, classification== 'Fetal Only')
d= select(d, chr, pos, EFF, REF, Beta_fetal, SE_fetal)
names(d)= c('CHR', 'POS', 'EFF', 'REF', 'beta', 'se')

d$beta= ifelse(d$REF> d$EFF, -1 * d$beta, d$beta)
d$ID= ifelse(d$REF> d$EFF, paste(d$CHR, d$POS, d$EFF, d$REF, sep= ':'),  paste(d$CHR, d$POS, d$REF, d$EFF, sep= ':'))

x= fread(snakemake@input[[2]], select= c('ID', 'BETA', 'SE'))

d= inner_join(d, x, by= 'ID')

inputMR= mr_input(bx = d$beta,   bxse= d$se,by = d$BETA, byse = d$SE)

z= mr_allmethods(inputMR)$Values
names(z)= c('method', 'estimate', 'se', 'lo95', 'up95', 'pvalue')
z= filter(z, method== 'IVW' | method== 'MR-Egger')
z$outcome= 'UL'

d= fread(snakemake@input[[3]])

moms= fread(snakemake@input[[4]])
names(moms)= c('IID', 'moms_GRS')
fets= fread(snakemake@input[[5]])
names(fets)= c('IID', 'fets_GRS')
dads= fread(snakemake@input[[6]])
names(dads)= c('IID', 'dads_GRS')

d= inner_join(d, moms, by= c('Mother'= 'IID')) %>% inner_join(., fets, by= c('Child'= 'IID')) %>% inner_join(., dads, by= c('Father'= 'IID'))

d= filter(d, spont== 1, !is.na(SVLEN_UL_DG), !is.na(SVLEN_SM_DG))

m1= lm(SVLEN_UL_DG~ moms_GRS + fets_GRS + dads_GRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + cohort + PARITY, d)

x= data.frame(summary(m1)$coefficients)[2:4, ]
names(x)= c('estimate', 'se', 'z', 'pvalue')
x$method= row.names(x)

ci= data.frame(confint(m1))[2:4, ]
names(ci)= c('lo95', 'up95')
ci$method= row.names(ci)
x= inner_join(x, ci, by= 'method')
x$outcome= 'UL'
x$sample_size= length(m1$resid)

m1= lm(SVLEN_SM_DG~ moms_GRS + fets_GRS + dads_GRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + cohort + PARITY, d)

x1= data.frame(summary(m1)$coefficients)[2:4, ]
names(x1)= c('estimate', 'se', 'z', 'pvalue')
x1$method= row.names(x1)

ci= data.frame(confint(m1))[2:4, ]
names(ci)= c('lo95', 'up95')
ci$method= row.names(ci)
x1= inner_join(x1, ci, by= 'method')
x1$outcome= 'LMP'
x1$sample_size= length(m1$resid)

x= rbind(x, x1)


d= fread(snakemake@input[[3]])

h1= fread(snakemake@input[[7]])
names(h1)= c('h1', 'PREG_ID')
h2= fread(snakemake@input[[8]])
names(h2)= c('h2', 'PREG_ID')
h3= fread(snakemake@input[[9]])
names(h3)= c('h3', 'PREG_ID')
h4= fread(snakemake@input[[10]])
names(h4)= c('h4', 'PREG_ID')

d= inner_join(h1, h2, by= 'PREG_ID') %>% inner_join(., h3, by= 'PREG_ID') %>% inner_join(., h4, by= 'PREG_ID') %>% inner_join(., d, by= 'PREG_ID')

d= filter(d, spont== 1, !is.na(SVLEN_UL_DG), !is.na(SVLEN_SM_DG))

m1= lm(SVLEN_UL_DG~ h1 + h2 + h3 + h4 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + cohort + PARITY, d)

h= data.frame(summary(m1)$coefficients)[2:4, ]
names(h)= c('estimate', 'se', 'z', 'pvalue')
h$method= row.names(h)

ci= data.frame(confint(m1))[2:4, ]
names(ci)= c('lo95', 'up95')
ci$method= row.names(ci)
h= inner_join(h, ci, by= 'method')
h$outcome= 'UL'
h$sample_size= length(m1$resid)

m1= lm(SVLEN_SM_DG~ h1 + h2 + h3 + h4 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + cohort + PARITY, d)

h1= data.frame(summary(m1)$coefficients)[2:4, ]
names(h1)= c('estimate', 'se', 'z', 'pvalue')
h1$method= row.names(h1)

ci= data.frame(confint(m1))[2:4, ]
names(ci)= c('lo95', 'up95')
ci$method= row.names(ci)
h1= inner_join(h1, ci, by= 'method')
h1$outcome= 'LMP'
h1$sample_size= length(m1$resid)

h= rbind(h, h1)


x= bind_rows(x, h)

x= bind_rows(z, x)

fwrite(x, snakemake@output[[1]], sep= '\t')

