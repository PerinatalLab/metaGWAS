library(dplyr)
library(data.table)
library(broom)


d= fread(snakemake@input[[1]])

d= filter(d, !duplicated(IID), fam_member == 'Mother')

m1= (lm(MnTGA_PGS ~ fetal_effect_PGS + maternal_effect_PGS, d))
cf1= data.frame(confint(m1))[2:3, ]
names(cf1)= c('lo95', 'up95')
cf1$exposure= rownames(cf1)

m1= tidy(m1)[2:3, ]
m1$outcome= 'Maternal non-transmitted PGS'

m1= inner_join(m1, cf1, by= c('term'= 'exposure'))

d= fread(snakemake@input[[2]])

d= filter(d, !duplicated(IID), fam_member == 'Mother')

m2= (lm(MTGA_PGS ~ fetal_effect_PGS + maternal_effect_PGS, d))

cf2= data.frame(confint(m2))[2:3, ]
names(cf2)= c('lo95', 'up95')
cf2$exposure= rownames(cf2)

m2= tidy(m2)[2:3, ]
m2$outcome= 'Maternal transmitted PGS'

m2= inner_join(m2, cf2, by= c('term'= 'exposure'))

d= fread(snakemake@input[[3]])

d= filter(d, !duplicated(IID), fam_member == 'Mother')

m3= (lm(PTGA_PGS ~ fetal_effect_PGS + maternal_effect_PGS, d))

cf3= data.frame(confint(m3))[2:3, ]
names(cf3)= c('lo95', 'up95')
cf3$exposure= rownames(cf1)

m3= tidy(m3)[2:3, ]
m2$outcome= 'Maternal transmitted PGS'

m3= inner_join(m3, cf3, by= c('term'= 'exposure'))
m3$outcome= 'Paternal transmitted PGS'


d= do.call('rbind', list(m1, m2, m3))

fwrite(d, snakemake@output[[1]], sep= '\t')

