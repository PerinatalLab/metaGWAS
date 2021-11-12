library(data.table)
library(tidyr)
library(dplyr)

top_pos= 123112292
dist_lim= 1.5*10**6


d= fread(snakemake@input[[1]], select= c('ID', 'CHR', 'POS', 'EAF', 'BETA', 'SE', 'TOTALSAMPLESIZE'))
d= filter(d, CHR== 3, POS>= top_pos - dist_lim,  POS<= top_pos + dist_lim)

d$MAF= ifelse(d$EAF> 0.5, 1 - d$EAF, d$EAF)
bw= fread(snakemake@input[[2]], select= c('ID', 'CHR', 'POS', 'EAF', 'BETA', 'SE', 'N'))
names(bw)= c('ID', 'CHR', 'POS', 'EAF', 'BETA', 'SE', 'TOTALSAMPLESIZE')
bw= filter(bw, CHR== 3, POS>= top_pos - dist_lim,  POS<= top_pos + dist_lim)

bw$MAF= ifelse(bw$EAF>0.5, 1 - bw$EAF, bw$EAF)

fin_map= fread(snakemake@input[[3]], select= c('ID', 'CHR', 'POS'))

names(fin_map)= c('ID_hg38', 'CHR', 'POS')

fin_map= filter(fin_map, CHR== 'chr3', POS>= top_pos - dist_lim,  POS<= top_pos + dist_lim)
fin_map= separate(fin_map, ID_hg38, into= c('chr_hg38', 'pos_hg38', 'ref', 'eff'), sep= ':', remove= FALSE)

pos_hg38= as.numeric(fin_map$pos_hg38)

fin_map$CHR= gsub('chr', '', fin_map$CHR)
fin_map$ID= paste(fin_map$CHR, fin_map$POS, fin_map$ref, fin_map$eff, sep= ':')

d= inner_join(d, fin_map[, c('ID', 'ID_hg38')], by= 'ID')
bw= inner_join(bw, fin_map[, c('ID', 'ID_hg38')], by= 'ID')


fwrite(d, snakemake@output[[1]], sep= '\t')
fwrite(bw, snakemake@output[[2]], sep= '\t')
