library(data.table)
library(dplyr)

q1= fread(snakemake@input[[1]])
q1= select(q1, PREG_ID_315, AA87, AA85, AA88)
names(q1)= c('PREG_ID', 'height', 'weight', 'father_height')

q1$BMI= with(q1, weight / ((height/ 100)**2))
q1= filter(q1, abs(scale(height))< 4, abs(scale(weight))< 4)
q1$PREG_ID= as.character(q1$PREG_ID)

fets= fread(snakemake@input[[2]])
moms= fread(snakemake@input[[3]])
dads= fread(snakemake@input[[4]])


format_haps= function(hap){
variants= paste(hap$chr, hap$pos, hap$ref, hap$eff, sep =':')
ids= names(hap)[5:ncol(hap)]
hap= as.data.frame(t(hap[, 5:ncol(hap)]))
names(hap)= variants
hap$IID= ids
return(hap)
}

fets= format_haps(fets)
moms= format_haps(moms)
dads= format_haps(dads)

names(fets)= c('BW_fets', 'GA_fets', 'Child')
names(moms)= c('BW_moms', 'GA_moms', 'Mother')
names(dads)= c('BW_dads', 'GA_dads', 'Father')

pheno= fread(snakemake@input[[5]])

pheno$PREG_ID= as.character(pheno$PREG_ID)


d= left_join(pheno, fets, by= 'Child') %>% left_join(., moms, by= 'Mother') %>% left_join(., dads, by= 'Father')
d= left_join(d, q1, by= 'PREG_ID')

m1= (lm(SVLEN_UL_DG~ GA_fets + GA_moms + GA_dads + BW_fets + BW_moms + BW_dads + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PARITY + cohort, filter(d, spont== 1)))


x= data.frame((summary(m1)$coefficients)[2:7, ])
names(x)= c('beta', 'se', 'z', 'pvalue')

x$exposure= row.names(x)

m1= (lm(SVLEN_UL_DG~ GA_fets + GA_moms + GA_dads + BW_fets + BW_moms + BW_dads + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PARITY + cohort + VEKT, filter(d, spont== 1)))


x1= data.frame((summary(m1)$coefficients)[2:7, ])
names(x1)= c('beta_bw', 'se_bw', 'z_bw', 'pvalue_bw')

x1$exposure= row.names(x1)

x= inner_join(x, x1, by= 'exposure')

x$locus= ifelse(grepl('BW', x$exposure), 'rs11708067 (birth weight)', 'rs28654158 (gestational duration)')
x$genome= ifelse(grepl('moms', x$exposure), 'Maternal', ifelse(grepl('fets', x$exposure), 'Fetal', 'Paternal'))


fwrite(x, snakemake@output[[1]], sep= '\t')



cat("
h1= fread('/mnt/work2/pol/metaGWAS/ADCY5/processed_data/haplotypes/h1_PREG_ID')
h2= fread('/mnt/work2/pol/metaGWAS/ADCY5/processed_data/haplotypes/h2_PREG_ID')
h3= fread('/mnt/work2/pol/metaGWAS/ADCY5/processed_data/haplotypes/h3_PREG_ID')
h4= fread('/mnt/work2/pol/metaGWAS/ADCY5/processed_data/haplotypes/h4_PREG_ID')

h1= format_haps(h1)
h2= format_haps(h2)
h3= format_haps(h3)
h4= format_haps(h4)

names(h1)= c('BW_h1', 'GA_h1', 'PREG_ID')
names(h2)= c('BW_h2', 'GA_h2', 'PREG_ID')
names(h3)= c('BW_h3', 'GA_h3', 'PREG_ID')
names(h4)= c('BW_h4', 'GA_h4', 'PREG_ID')

x= left_join(pheno, h1, by= 'PREG_ID') %>% left_join(., h2, by= 'PREG_ID')%>% left_join(., h3, by= 'PREG_ID') %>% left_join(., h4, by= 'PREG_ID')

x= left_join(x, q1, by= 'PREG_ID')


m1= (lm(SVLEN_UL_DG~ GA_fets + GA_moms + GA_dads + BW_fets + BW_moms + BW_dads + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PARITY + cohort, filter(d, spont== 1)))


x= data.frame((summary(m1)$coefficients)[2:7, ])
names(x)= c('beta', 'se', 'z', 'pvalue')

x$exposure= row.names(x)

m1= (lm(SVLEN_UL_DG~ GA_fets + GA_moms + GA_dads + BW_fets + BW_moms + BW_dads + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PARITY + cohort + VEKT, filter(d, spont== 1)))


x1= data.frame((summary(m1)$coefficients)[2:7, ])
names(x1)= c('beta_bw', 'se_bw', 'z_bw', 'pvalue_bw')

x1$exposure= row.names(x1)

x= inner_join(x, x1, by= 'exposure')

x$locus= ifelse(grepl('BW', x$exposure), 'rs11708067 (birth weight)', 'rs28654158 (gestational duration)')
x$genome= ifelse(grepl('moms', x$exposure), 'Maternal', ifelse(grepl('fets', x$exposure), 'Fetal', 'Paternal'))


fwrite(x, snakemake@output[[1]], sep= '\t')
")
