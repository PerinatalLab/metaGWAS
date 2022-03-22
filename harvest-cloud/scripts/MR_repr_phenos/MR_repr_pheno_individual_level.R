library(dplyr)
library(data.table)


d= fread(snakemake@input[[1]])

m1= lm(SVLEN_UL_DG~ h1 + h2 + h3  + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + cohort + PARITY, d)

h= data.frame(summary(m1)$coefficients)[2:4, ]
names(h)= c('estimate', 'se', 'z', 'pvalue')
h$method= row.names(h)

ci= data.frame(confint(m1))[2:4, ]
names(ci)= c('lo95', 'up95')
ci$method= row.names(ci)
h= inner_join(h, ci, by= 'method')
h$outcome= 'UL'

h$exposure= snakemake@wildcards[['repr_trait']]
h$n= length(m1$resid)

m1= lm(SVLEN_UL_DG~ h1 + h2 + h3 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + cohort + PARITY, filter(d, spont== 1))

x= data.frame(summary(m1)$coefficients)[2:4, ]
names(x)= c('estimate', 'se', 'z', 'pvalue')
x$method= row.names(x)

ci= data.frame(confint(m1))[2:4, ]
names(ci)= c('lo95', 'up95')
ci$method= row.names(ci)
x= inner_join(x, ci, by= 'method')
x$outcome= 'UL_spont'

x$exposure= snakemake@wildcards[['repr_trait']]
x$n= length(m1$resid)

h= rbind(h, x)

fwrite(h, snakemake@output[[1]], sep= '\t')

