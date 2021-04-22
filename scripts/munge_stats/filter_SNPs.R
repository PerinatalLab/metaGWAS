library(data.table)
library(dplyr)

d= fread(snakemake@input[[1]])

x1= nrow(d)

d= arrange(d, CHR, POS, EFF, REF)

hrc= fread(snakemake@input[[2]], header=T)

d= inner_join(d, hrc, by= 'ID')
rm(hrc)
d$EAF= ifelse(is.na(d$EAF), d$eaf, d$EAF)

d$BETA= ifelse(d$REF> d$EFF, -1 * d$BETA, d$BETA)
d$EAF= ifelse(d$REF> d$EFF, 1 - d$EAF, d$EAF)

d[d$REF>d$EFF, c("REF", "EFF")]= d[d$REF > d$EFF, c("EFF", "REF")]

d$MAF= ifelse(d$EAF>0.5, 1- d$EAF, d$EAF)

d= filter(d, pvalue>0, pvalue<1, MAF>0.005, SE>0)

d= filter(d, (MAF * 2 * N) > 6)

d$maf= ifelse(d$eaf> 0.5, 1 - d$eaf, d$eaf)

d= filter(d, abs(MAF - maf) < 0.2)

d= select(d, -c(MAF, ID, eaf))

x2= nrow(d)

write.table(d, snakemake@output[[1]], col.names= T, row.names=F, sep= '\t', quote= F)

cohort= gsub('_temp.txt', '', unlist(strsplit(snakemake@input[[1]], '/'))[[9]])

cat(c(cohort, '\t', x1, '\t', x2, '\n'), file= snakemake@output[[2]])
