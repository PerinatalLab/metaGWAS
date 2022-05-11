library(dplyr)
library(data.table)

format_haps= function(hap){
variants= paste(hap$chr, hap$pos, hap$ref, hap$eff, sep =':')
ids= names(hap)[5:ncol(hap)]
hap= as.data.frame(t(hap[, 5:ncol(hap)]))
names(hap)= variants
hap$PREG_ID= as.numeric(ids)
return(hap)
}

write(paste('snp', 'n', 'freq_h1', 'freq_h2', 'freq_h3', 'beta_h1', 'se_h1', 'pvalue_h1', 'beta_h2', 'se_h2', 'pvalue_h2', 'beta_h3', 'se_h3', 'pvalue_h3', sep= '\t'), snakemake@output[[1]], append= T)


h1= fread(snakemake@input[[1]], h= T)
h2= fread(snakemake@input[[2]], h= T)
h3= fread(snakemake@input[[3]], h= T)

h1= format_haps(h1)
names(h1)= c('h1', 'PREG_ID')
h2= format_haps(h2)
names(h2)= c('h2', 'PREG_ID')
h3= format_haps(h3)
names(h3)= c('h3', 'PREG_ID')

pheno= fread(snakemake@input[[5]])
pheno= filter(pheno, spont== 1)

h= inner_join(h1, h2, by= 'PREG_ID') %>% inner_join(., h3, by= 'PREG_ID')
pheno$PREG_ID= as.numeric(pheno$PREG_ID)

d= inner_join(h, pheno, by= 'PREG_ID')

m1= lm(SVLEN_UL_DG~ h1 + h2 + h3 + PARITY + cohort + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, d)

n= length(resid(m1))
h1= fread(snakemake@input[[1]], h= T)

snp= paste(h1$chr, h1$pos, h1$ref, h1$eff, sep =':')
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

freq_h1= mean(d$h1, na.rm= T)
freq_h2= mean(d$h2, na.rm= T)
freq_h3= mean(d$h3, na.rm= T)

results= paste(snp, n, freq_h1, freq_h2, freq_h3, beta_h1, se_h1, pvalue_h1, beta_h2, se_h2, pvalue_h2, beta_h3, se_h3, pvalue_h3, sep= '\t')
write(results, file= snakemake@output[[1]], append=TRUE)


