library("dplyr")
library("knitr")
library("tidyr")
library(cowplot)
library(ggrepel)
library("data.table")
library('showtext')



colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


font_add("arial", "arial.ttf", bold= 'arial_bold.ttf')

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)

d= fread(snakemake@input[[1]])

d= filter(d, MarkerName!= '6:32595083:G:T')

d= separate(d, MarkerName, into= c('CHR', 'POS', 'REF', 'EFF'), sep= ':')
d$beta_h1= with(d, ifelse(Allele2 > Allele1, -1 * beta_h1, beta_h1))
d$beta_h2= with(d, ifelse(Allele2 > Allele1, -1 * beta_h2, beta_h2))
d$beta_h3= with(d, ifelse(Allele2 > Allele1, -1 * beta_h3, beta_h3))


d$ID= with(d, ifelse(REF> EFF, paste(CHR, POS, EFF, REF, sep= ':'), paste(CHR, POS, REF, EFF, sep= ':')))


top= fread(snakemake@input[[2]])
top= filter(top, ID!= '6:32595083:G:T')
ids= pull(top, ID)
ids= c('3:156697097:A:G', '5:158058432:G:T', ids)

d= filter(d, ID %in% ids)

x= fread(snakemake@input[[3]], select= c('ID', 'BETA'))

x= filter(x, ID %in% ids)

d= inner_join(d, x, by= 'ID')

d$beta_h2= with(d, ifelse(BETA< 0, -1 * beta_h2,  beta_h2))
d$beta_h3= with(d, ifelse(BETA< 0, -1 * beta_h3, beta_h3))
d$beta_h1= with(d, ifelse(BETA< 0, -1 * beta_h1, beta_h1))
d$BETA= with(d, ifelse(BETA<0, -1 * BETA, BETA))

d$lo95_h1= d$beta_h1 - 1.96 * d$se_h1
d$up95_h1= d$beta_h1 + 1.96 * d$se_h1

d$lo95_h2= d$beta_h2 - 1.96 * d$se_h2
d$up95_h2= d$beta_h2 + 1.96 * d$se_h2

d$lo95_h3= d$beta_h3 - 1.96 * d$se_h3
d$up95_h3= d$beta_h3 + 1.96 * d$se_h3

d$sig_origin= with(d, ifelse(pvalue_h1<0.05 & pvalue_h2>= 0.05 & pvalue_h3 >= 0.05, 'Maternal and/or fetal', ifelse(pvalue_h1<0.05 & pvalue_h2<0.05 & pvalue_h3>= 0.05, 'Maternal', ifelse(pvalue_h1>= 0.05 & pvalue_h2< 0.05 & pvalue_h3>= 0.05, 'Maternal', ifelse(pvalue_h1>= 0.05 & pvalue_h2< 0.05 & pvalue_h3< 0.05, 'Maternal and fetal', 'Unclassified')))))

d$sig_origin= factor(d$sig_origin, levels= c('Maternal', 'Fetal', 'Maternal and fetal', 'Maternal and/or fetal', 'Unclassified'))

d$sig_h2= with(d, ifelse(pvalue_h2<0.05, '2', '1'))
d$sig_h1= with(d, ifelse(pvalue_h1<0.05, '2', '1'))
d$sig_h3= with(d, ifelse(pvalue_h3<0.05, '2', '1'))

p1= ggplot(d, aes(beta_h2, BETA, colour= sig_h2)) +
geom_point(size= 0.5) +
#geom_errorbarh(data= filter(d, (lo95_h2 >0 & up95_h2>0) | (lo95_h2<0 & up95_h2 <0)), aes(xmax = lo95_h2, xmin = up95_h2), size= 0.05) +
theme_cowplot(font_size= 8) +
scale_colour_manual(values= c('grey', colorBlindBlack8[8]), guide= F) +
geom_vline(xintercept= 0, colour= colorBlindBlack8[1], linetype= 'dashed', size= 0.2, alpha= 0.6) +
geom_hline(yintercept= 0, colour= colorBlindBlack8[1], linetype= 'dashed', size= 0.2, alpha= 0.6) +
xlab('Effect size maternal \nnon-transmitted alleles, days') +
ylab('Effect size maternal genome, days') 
#scale_x_continuous(breaks = round(seq(-1.5, 3, by= 0.5), 1)) +
#  scale_y_continuous(breaks = round(seq(-1.5, 3, by= 0.5), 1))


ggsave(snakemake@output[[1]], plot= p1, width= 60, height= 60, units= 'mm', dpi= 300)

print('plot1')
p1= ggplot(d, aes(beta_h3, BETA, colour= sig_h3)) +
geom_point(size= 0.5) +
#geom_errorbarh(data= filter(d, (lo95_h3 >0 & up95_h3>0) | (lo95_h3<0 & up95_h3 <0)), aes(xmax = lo95_h3, xmin = up95_h3), size= 0.05) +
theme_cowplot(font_size= 8) +
scale_colour_manual(values= c('grey', colorBlindBlack8[8]), guide= F) +
geom_vline(xintercept= 0, colour= colorBlindBlack8[1], linetype= 'dashed', size= 0.2, alpha= 0.6) +
geom_hline(yintercept= 0, colour= colorBlindBlack8[1], linetype= 'dashed', size= 0.2, alpha= 0.6) +
xlab('Effect size paternal \ntransmitted alleles, days') +
ylab('Effect size maternal genome, days') 
#scale_x_continuous(breaks = round(seq(-1.5, 3, by= 0.5), 1)) +
#  scale_y_continuous(breaks = round(seq(-1.5, 3, by= 0.5), 1))

ggsave(snakemake@output[[2]], plot= p1, width= 60, height= 60, units= 'mm', dpi= 300)

print('plot2')
p1= ggplot(d, aes(beta_h1, BETA, colour= sig_h1)) +
geom_point(size= 0.5) +
#geom_errorbarh(data= filter(d, (lo95_h3 >0 & up95_h3>0) | (lo95_h3<0 & up95_h3 <0)), aes(xmax = lo95_h3, xmin = up95_h3), size= 0.05) +
theme_cowplot(font_size= 8) +
scale_colour_manual(values= c('grey', colorBlindBlack8[8]), guide= F) +
geom_vline(xintercept= 0, colour= colorBlindBlack8[1], linetype= 'dashed', size= 0.2, alpha= 0.6) +
geom_hline(yintercept= 0, colour= colorBlindBlack8[1], linetype= 'dashed', size= 0.2, alpha= 0.6) +
xlab('Effect size maternal \ntransmitted alleles, days') +
ylab('Effect size maternal genome, days') 
#scale_x_continuous(breaks = round(seq(-1.5, 3, by= 0.5), 1)) +
#  scale_y_continuous(breaks = round(seq(-1.5, 3, by= 0.5), 1))

ggsave(snakemake@output[[3]], plot= p1, width= 60, height= 60, units= 'mm', dpi= 300)


fwrite(d, snakemake@output[[4]], sep= '\t')
