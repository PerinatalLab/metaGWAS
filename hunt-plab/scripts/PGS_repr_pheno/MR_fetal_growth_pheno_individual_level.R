library(data.table)
library(dplyr)


MnT= fread(snakemake@input[[1]])

MT= fread(snakemake@input[[2]])

PT= fread(snakemake@input[[3]])

pheno= fread(snakemake@input[[4]])
pheno$PREG_ID= paste('PREG', pheno$MOR_PID, pheno$BARN_PID, pheno$FAR_PID, sep= '_')
trio= fread(snakemake@input[[5]])

pheno= filter(pheno, PREG_ID %in% trio$PREG_ID)
print(nrow(pheno))
d= inner_join(MnT, MT, by= 'PREG_ID') %>% inner_join(., PT, by= 'PREG_ID') %>% inner_join(., pheno, by= 'PREG_ID')

d= arrange(d, desc(FAAR))
d= filter(d, !duplicated(MOR_PID))
d= filter(d, !duplicated(FAR_PID))

m1= lm(SVLEN_SM_DG~ MT + MnT + PT + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PARITY0 + FAAR, d)

x= data.frame(summary(m1)$coefficients)[2:4, ]
names(x)= c('estimate', 'se', 'z', 'pvalue')
x$method= row.names(x)

ci= data.frame(confint(m1))[2:4, ]
names(ci)= c('lo95', 'up95')
ci$method= row.names(ci)
x= inner_join(x, ci, by= 'method')
x$exposure= snakemake@wildcards[['repr_trait']]

x$n= length(m1$resid)





fwrite(x, snakemake@output[[1]], sep= '\t')
