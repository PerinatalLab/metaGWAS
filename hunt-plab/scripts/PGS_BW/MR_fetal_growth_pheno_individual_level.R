library(data.table)
library(dplyr)


MnT= fread(snakemake@input[[1]])
names(MnT)= c('h2', 'PREG_ID')

MT= fread(snakemake@input[[2]])
names(MT)= c('h1', 'PREG_ID')

PT= fread(snakemake@input[[3]])
names(PT)= c('h3', 'PREG_ID')

pheno= fread(snakemake@input[[4]])
pheno$PREG_ID= paste('PREG', pheno$MOR_PID, pheno$BARN_PID, pheno$FAR_PID, sep= '_')
trio= fread(snakemake@input[[5]])

pheno= filter(pheno, PREG_ID %in% trio$PREG_ID)
print(nrow(pheno))
d= inner_join(MnT, MT, by= 'PREG_ID') %>% inner_join(., PT, by= 'PREG_ID') %>% inner_join(., pheno, by= 'PREG_ID')

d= filter(d, !duplicated(MOR_PID))
d= filter(d, !duplicated(FAR_PID))

m1= lm(SVLEN_SM_DG~ h1 + h2 + h3 +PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PARITY0 + FAAR, d)

x= data.frame(summary(m1)$coefficients)[2:4, ]
names(x)= c('estimate', 'se', 'z', 'pvalue')
x$method= row.names(x)

ci= data.frame(confint(m1))[2:4, ]
names(ci)= c('lo95', 'up95')
ci$method= row.names(ci)
x= inner_join(x, ci, by= 'method')
x$outcome= 'SM_spont'

x$exposure= 'fetal_growth'
x$n= length(m1$resid)





fwrite(x, snakemake@output[[1]], sep= '\t')
