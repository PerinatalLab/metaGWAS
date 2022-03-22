library(dplyr)
library(data.table)


d= fread(snakemake@input[[1]])

m1= lm(SVLEN_UL_DG~ h1 + h2 + h3 + h4 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + cohort + PARITY, d)

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

m1= lm(SVLEN_UL_DG~ h1 + h2 + h3 + h4 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + cohort + PARITY, filter(d, spont== 1))

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

z= rbind(h, x)

d= fread(snakemake@input[[1]])

m1= lm(SVLEN_SM_DG~ h1 + h2 + h3 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + cohort + PARITY, d)

h= data.frame(summary(m1)$coefficients)[2:4, ]
names(h)= c('estimate', 'se', 'z', 'pvalue')
h$method= row.names(h)

ci= data.frame(confint(m1))[2:4, ]
names(ci)= c('lo95', 'up95')
ci$method= row.names(ci)
h= inner_join(h, ci, by= 'method')
h$outcome= 'SM'

h$exposure= snakemake@wildcards[['repr_trait']]
h$n= length(m1$resid)

print(cor(d$h3, d$h4, use= 'complete'))
print(cor(d$h2, d$h1, use= 'complete'))
print(cor(d$h2, d$h4, use= 'complete'))

m1= lm(SVLEN_SM_DG~ h1 + h2 + h3 + h4 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + cohort + PARITY, filter(d, spont== 1))

x= data.frame(summary(m1)$coefficients)[2:4, ]
names(x)= c('estimate', 'se', 'z', 'pvalue')
x$method= row.names(x)

ci= data.frame(confint(m1))[2:4, ]
names(ci)= c('lo95', 'up95')
ci$method= row.names(ci)
x= inner_join(x, ci, by= 'method')
x$outcome= 'SM_spont'

x$exposure= snakemake@wildcards[['repr_trait']]
x$n= length(m1$resid)
z1= rbind(h, x)

df= rbind(z, z1)

fwrite(df, snakemake@output[[1]], sep= '\t')

