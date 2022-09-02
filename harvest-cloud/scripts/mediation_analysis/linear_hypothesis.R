library(data.table)
library(dplyr)
library(tidyr)
library(car)

format_haps= function(hap){
print(dim(hap))
variants= paste(hap$chr, hap$pos, hap$ref, hap$eff, sep =':')
ids= names(hap)[5:ncol(hap)]
hap= as.data.frame(t(hap[, 5:ncol(hap)]))
names(hap)= variants
hap$PREG_ID= as.character(as.integer(ids))
return(hap)
}

h1= fread(snakemake@input[[1]])
h2= fread(snakemake@input[[2]])
h3= fread(snakemake@input[[3]])
h4= fread(snakemake@input[[4]])

h1= format_haps(h1)
h2= format_haps(h2)
h3= format_haps(h3)
h4= format_haps(h4)

pheno= fread(snakemake@input[[5]])
#pheno= filter(pheno, spont== 1)

pheno$PREG_ID= as.character(as.integer(pheno$PREG_ID))

print(nrow(pheno))
write(paste('snp', 'n', 'freq_h1', 'freq_h2', 'freq_h3', 'beta_h1', 'se_h1', 'pvalue_h1', 'beta_h2', 'se_h2', 'pvalue_h2', 'beta_h3', 'se_h3', 'pvalue_h3', 'effect', 'beta_h1_GA', 'se_h1_GA', 'pvalue_h1_GA', 'beta_h2_GA', 'se_h2_GA', 'pvalue_h2_GA', 'beta_h3_GA', 'se_h3_GA', 'pvalue_h3_GA', sep= '\t'), snakemake@output[[1]], append= T)

results_list= lapply(names(h1)[1:(length(names(h1))-1)], function(snp) {

if (grepl('X', snp)){ 
print('Not sure how to handle chromosome X.')

} else {

h1_temp= h1[, c('PREG_ID', snp)]
h2_temp= h2[, c('PREG_ID', snp)]
h3_temp= h3[, c('PREG_ID', snp)]
h4_temp= h4[, c('PREG_ID', snp)]

names(h1_temp)= c('PREG_ID', 'h1')
names(h2_temp)= c('PREG_ID', 'h2')
names(h3_temp)= c('PREG_ID', 'h3')
names(h4_temp)= c('PREG_ID', 'h4')

d= inner_join(pheno, h1_temp, by= 'PREG_ID') %>% inner_join(., h2_temp, by= 'PREG_ID') %>% inner_join(., h3_temp, by= 'PREG_ID')

d= filter(d, !is.na(VEKT), !is.na(SVLEN_UL_DG))

m1= lm(VEKT~ h1 + h2 + h3 + cohort + PC1 + PC2 + PC3 + PC4 + PC5 + PC6, d)

n= length(resid(m1))
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

m2= lm(VEKT~ h1 + h2 + h3 + cohort + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + SVLEN_UL_DG, d)

coefs= summary(m2)$coefficients[2:5,]
beta_h12= coefs[1,1]
se_h12= coefs[1,2]
pvalue_h12= coefs[1,4]
beta_h22= coefs[2,1]
se_h22= coefs[2,2]
pvalue_h22= coefs[2,4]
beta_h32= coefs[3,1]
se_h32= coefs[3,2]
pvalue_h32= coefs[3,4]


results= paste(snp, n, freq_h1, freq_h2, freq_h3, beta_h1, se_h1, pvalue_h1, beta_h2, se_h2, pvalue_h2, beta_h3, se_h3, pvalue_h3, snakemake@wildcards[["effect_origin_BW"]], beta_h12, se_h12, pvalue_h12, beta_h22, se_h22, pvalue_h22, beta_h32, se_h32, pvalue_h32, sep= '\t')

write(results, file= snakemake@output[[1]], append=TRUE)

}
}
)



