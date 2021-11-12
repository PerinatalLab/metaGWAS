library(data.table)
library(dplyr)


hrc= fread(snakemake@input[[1]])


funk= function(infile){
d= fread(infile)

print(paste('Filtering file: ', infile))

d= arrange(d, CHR, POS, EFF, REF)

d= filter(d, pvalue< 1, pvalue>0)

d$pval= pnorm(-abs(d$BETA / d$SE)) * 2

d= filter(d, (abs(-log10(pvalue) - -log10(pval)) / -log10(pval)) * 100 <= 10)

d$ID= with(d, ifelse(REF> EFF, paste(CHR, POS, EFF, REF, sep= ':'), paste(CHR, POS, REF, EFF, sep= ':')))

d$SNP= with(d, ifelse(grepl('I', ID), paste(ID, 'INDEL', sep= ':'), paste(ID, 'SNP', sep= ':')))


print(str(d))

d= inner_join(d, hrc, by= 'ID')
d$EAF= ifelse(is.na(d$EAF), d$eaf, d$EAF)

d$BETA= ifelse(d$REF> d$EFF, -1 * d$BETA, d$BETA)
d$EAF= ifelse(d$REF> d$EFF, 1 - d$EAF, d$EAF)

d[d$REF>d$EFF, c("REF", "EFF")]= d[d$REF > d$EFF, c("EFF", "REF")]

d$MAF= ifelse(d$EAF>0.5, 1- d$EAF, d$EAF)
d= filter(d, MAF>= 0.005)

d= filter(d, pvalue>0, pvalue<1, MAF>=0.005, SE>0)

d= filter(d, (MAF * 2 * N) > 6)

d$maf= ifelse(d$eaf> 0.5, 1 - d$eaf, d$eaf)

d= filter(d, abs(maf - MAF) < 0.2)

if (grepl('GAraw/Viva', infile)){



d$EAF= with(d, ifelse(abs(eaf - EAF)> 0.2, 1 - EAF, EAF))
d$BETA= with(d, ifelse(abs(eaf - EAF)> 0.2, -1 * BETA, BETA))
}


if (grepl('GAnrm/Viva', infile)){



d$EAF= with(d, ifelse(abs(eaf - EAF)> 0.2, 1 - EAF, EAF))
d$BETA= with(d, ifelse(abs(eaf - EAF)> 0.2, -1 * BETA, BETA))
}

if (grepl('postTerm/HUNT', infile)){



d$EAF= with(d, ifelse(abs(eaf - EAF)> 0.2, 1 - EAF, EAF))
d$BETA= with(d, ifelse(abs(eaf - EAF)> 0.2, -1 * BETA, BETA))
}


d= arrange(d, pvalue)
d= filter(d, !duplicated(ID))


d= select(d, -c(MAF, ID, eaf, pval))

x2= nrow(d)

d$STRAND= '+'

#outfile= paste0(snakemake@params[[1]], gsub('_temp.txt', '', unlist(strsplit(infile, '/'))[9]), '.txt')


fwrite(d, snakemake@output[[1]], sep= '\t')
}


#input_files= snakemake@input[grepl('sumstats', snakemake@input)]

lapply(snakemake@input[[2]], funk)

