library(data.table)
library(dplyr)
library(tidyr)
library(car)

format_haps= function(hap){
variants= paste(hap$chr, hap$pos, hap$ref, hap$eff, sep =':')
ids= names(hap)[5:ncol(hap)]
hap= as.data.frame(t(hap[, 5:ncol(hap)]))
names(hap)= variants
hap$PREG_ID= ids
return(hap)
}

h1= fread(snakemake@input[[1]])
h2= fread(snakemake@input[[2]])
h3= fread(snakemake@input[[3]])
#h4= fread(snakemake@input[[4]])

h1= format_haps(h1)
h2= format_haps(h2)
h3= format_haps(h3)
#h4= format_haps(h4)

pheno= fread(snakemake@input[[5]])
pheno= filter(pheno, spont== 1)

write( paste('snp', 'n', 'freq_h1', 'freq_h2', 'freq_h3', 'beta_h1', 'se_h1', 'pvalue_h1', 'beta_h2', 'se_h2', 'pvalue_h2', 'beta_h3', 'se_h3', 'pvalue_h3', 'pval_maternal', 'pval_fetal', 'pval_poe', sep= '\t'), snakemake@output[[1]], append= T)

results_list= lapply(names(h1)[1:(length(names(h1))-1)], function(snp) {
h1_temp= h1[, c('PREG_ID', snp)]
h2_temp= h2[, c('PREG_ID', snp)]
h3_temp= h3[, c('PREG_ID', snp)]
#h4_temp= h4[, c('PREG_ID', snp)]

names(h1_temp)= c('PREG_ID', 'h1')
names(h2_temp)= c('PREG_ID', 'h2')
names(h3_temp)= c('PREG_ID', 'h3')
#names(h4_temp)= c('PREG_ID', 'h4')

d= inner_join(pheno, h1_temp, by= 'PREG_ID') %>% inner_join(., h2_temp, by= 'PREG_ID') %>% inner_join(., h3_temp, by= 'PREG_ID')# %>% inner_join(., h4_temp, by= 'PREG_ID')

n= nrow(d)
#m1= lm(SVLEN_UL_DG~ h1 + h2 + h3 + h4 + PARITY0 + cohort, d)

m1= lm(SVLEN_UL_DG~ h1 + h2 + h3 + PARITY0 + cohort + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, d)
coefs= summary(m1)$coefficients[2:5,]
beta_h1= coefs[1,1]
se_h1= coefs[1,2]
pvalue_h1= coefs[1,4]
beta_h2= coefs[2,1]
se_h2= coefs[2,2]
pvalue_h2= coefs[2,4]
beta_h3= coefs[3,1]
se_h3= coefs[3,2]
pvalue_h3= coefs[3,4]
#beta_h4= coefs[4,1]
#se_h4= coefs[4,2]
#pvalue_h4= coefs[4,4]

freq_h1= mean(d$h1, na.rm= T)
freq_h2= mean(d$h2, na.rm= T)
freq_h3= mean(d$h3, na.rm= T)

pval_maternal= tryCatch(linearHypothesis(m1, 'h1 + h2 = h3')[['Pr(>F)']][2], warning= function(w){NA}, error= function(w) {NA})
pval_fetal= tryCatch(linearHypothesis(m1, 'h1 + h3 = h2')[['Pr(>F)']][2], warning= function(w){NA}, error= function(w) {NA})
pval_poe= tryCatch(linearHypothesis(m1, 'h1 - h2 = h3')[['Pr(>F)']][2], warning= function(w){NA}, error= function(w) {NA})

results= paste(snp, n, freq_h1, freq_h2, freq_h3, beta_h1, se_h1, pvalue_h1, beta_h2, se_h2, pvalue_h2, beta_h3, se_h3, pvalue_h3, pval_maternal, pval_fetal, pval_poe, sep= '\t')
write(results, file= snakemake@output[[1]], append=TRUE)
}
)



