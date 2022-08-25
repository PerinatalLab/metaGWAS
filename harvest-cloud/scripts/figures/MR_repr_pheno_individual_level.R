library(dplyr)
library(data.table)


d= fread(snakemake@input[[1]])

d= fread(snakemake@input[[1]])

h1= fread(snakemake@input[[5]])
names(h1)= c('h1', 'PREG_ID')
h2= fread(snakemake@input[[6]])
names(h2)= c('h2', 'PREG_ID')
h3= fread(snakemake@input[[7]])
names(h3)= c('h3', 'PREG_ID')
h4= fread(snakemake@input[[8]])
names(h4)= c('h4', 'PREG_ID')

d= inner_join(h1, h2, by= 'PREG_ID') %>% inner_join(., h3, by= 'PREG_ID') %>% inner_join(., h4, by= 'PREG_ID') %>% inner_join(., d, by= 'PREG_ID')

m1= lm(SVLEN_UL_DG~ h1 + h2 + h3 + h4 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + cohort + PARITY, d)

h= data.frame(summary(m1)$coefficients)[2:4, ]
names(h)= c('estimate', 'se', 'z', 'pvalue')
h$method= row.names(h)

ci= data.frame(confint(m1))[2:4, ]
names(ci)= c('lo95', 'up95')
ci$method= row.names(ci)
h= inner_join(h, ci, by= 'method')
h$outcome= 'UL'

m1= lm(SVLEN_SM_DG~ h1 + h2 + h3 + h4 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + cohort + PARITY, d)

h1= data.frame(summary(m1)$coefficients)[2:4, ]
names(h1)= c('estimate', 'se', 'z', 'pvalue')
h1$method= row.names(h1)

ci= data.frame(confint(m1))[2:4, ]
names(ci)= c('lo95', 'up95')
ci$method= row.names(ci)
h1= inner_join(h1, ci, by= 'method')
h1$outcome= 'LMP'
h= bind_rows(h, h1)


x= bind_rows(x, h)


x$exposure= unlist(strsplit(unlist(strsplit(snakemake@input[[2]], '/'))[9], '_GRS'))[1]

fwrite(x, snakemake@output[[1]], sep= '\t')

