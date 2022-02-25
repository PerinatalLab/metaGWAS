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

x= fread(snakemake@input[[2]], select= c('RSID', 'BETA'))


d= inner_join(d, x, by= c('rsid' = 'RSID'))

d$beta_MNT= with(d, ifelse(BETA< 0, -1 * beta_MNT,  beta_MNT))
d$beta_PT= with(d, ifelse(BETA< 0, -1 * beta_PT, beta_PT))
d$beta_MT= with(d, ifelse(BETA< 0, -1 * beta_MT, beta_MT))
d$BETA= with(d, ifelse(BETA<0, -1 * BETA, BETA))

d$lo95_MT= d$beta_MT - 1.96 * d$se_MT
d$up95_MT= d$beta_MT + 1.96 * d$se_MT

d$lo95_MNT= d$beta_MNT - 1.96 * d$se_MNT
d$up95_MNT= d$beta_MNT + 1.96 * d$se_MNT

d$lo95_PT= d$beta_PT - 1.96 * d$se_PT
d$up95_PT= d$beta_PT + 1.96 * d$se_PT

d$class_name= with(d, ifelse(class_name== 'MF SD', 'Maternal and fetal (same direction)', ifelse(class_name== 'Fetal MatT', 'Fetal effect, maternal transmitted only', ifelse(class_name== 'Maternal', 'Maternal', ifelse(class_name== 'Fetal', 'Fetal', ifelse(class_name== 'MF OD', 'Maternal and fetal (opposite direction)', ''))))))

p1= ggplot(d, aes(beta_MNT, BETA, colour= class_name)) +
geom_point(size= 0.5) +
#geom_errorbarh(data= filter(d, (lo95_h2 >0 & up95_h2>0) | (lo95_h2<0 & up95_h2 <0)), aes(xmax = lo95_h2, xmin = up95_h2), size= 0.05) +
theme_cowplot(font_size= 8) +
scale_colour_manual(values= c('grey', colorBlindBlack8[c(8, 2, 4, 3)])) +
geom_vline(xintercept= 0, colour= colorBlindBlack8[1], linetype= 'dashed', size= 0.2, alpha= 0.6) +
geom_hline(yintercept= 0, colour= colorBlindBlack8[1], linetype= 'dashed', size= 0.2, alpha= 0.6) +
xlab('Effect size maternal \nnon-transmitted alleles, days') +
ylab('Effect size maternal genome, days')
#theme(legend.direction = "horizontal", legend.position = "bottom")
#scale_x_continuous(breaks = round(seq(-1.5, 3, by= 0.5), 1)) +
#  scale_y_continuous(breaks = round(seq(-1.5, 3, by= 0.5), 1))


ggsave(snakemake@output[[1]], plot= p1, width= 60, height= 60, units= 'mm', dpi= 300)

print('plot1')
p1= ggplot(d, aes(beta_PT, BETA, colour= class_name)) +
geom_point(size= 0.5) +
#geom_errorbarh(data= filter(d, (lo95_h3 >0 & up95_h3>0) | (lo95_h3<0 & up95_h3 <0)), aes(xmax = lo95_h3, xmin = up95_h3), size= 0.05) +
theme_cowplot(font_size= 8) +
scale_colour_manual(values= c('grey', colorBlindBlack8[c(8, 2, 4, 3)])) +
geom_vline(xintercept= 0, colour= colorBlindBlack8[1], linetype= 'dashed', size= 0.2, alpha= 0.6) +
geom_hline(yintercept= 0, colour= colorBlindBlack8[1], linetype= 'dashed', size= 0.2, alpha= 0.6) +
xlab('Effect size paternal \ntransmitted alleles, days') +
ylab('Effect size maternal genome, days') 
#scale_x_continuous(breaks = round(seq(-1.5, 3, by= 0.5), 1)) +
#  scale_y_continuous(breaks = round(seq(-1.5, 3, by= 0.5), 1))

ggsave(snakemake@output[[2]], plot= p1, width= 60, height= 60, units= 'mm', dpi= 300)

print('plot2')
p1= ggplot(d, aes(beta_MT, BETA, colour= class_name)) +
geom_point(size= 0.5) +
#geom_errorbarh(data= filter(d, (lo95_h3 >0 & up95_h3>0) | (lo95_h3<0 & up95_h3 <0)), aes(xmax = lo95_h3, xmin = up95_h3), size= 0.05) +
theme_cowplot(font_size= 8) +
scale_colour_manual(values= c('grey', colorBlindBlack8[c(8, 2, 4, 3)]), guide= F) +
geom_vline(xintercept= 0, colour= colorBlindBlack8[1], linetype= 'dashed', size= 0.2, alpha= 0.6) +
geom_hline(yintercept= 0, colour= colorBlindBlack8[1], linetype= 'dashed', size= 0.2, alpha= 0.6) +
xlab('Effect size maternal \ntransmitted alleles, days') +
ylab('Effect size maternal genome, days')
#scale_x_continuous(breaks = round(seq(-1.5, 3, by= 0.5), 1)) +
#  scale_y_continuous(breaks = round(seq(-1.5, 3, by= 0.5), 1))

ggsave(snakemake@output[[3]], plot= p1, width= 60, height= 60, units= 'mm', dpi= 300)

p1= ggplot(d, aes(beta_MNT, BETA, colour= class_name)) +
geom_point(size= 0.5) +
#geom_errorbarh(data= filter(d, (lo95_h2 >0 & up95_h2>0) | (lo95_h2<0 & up95_h2 <0)), aes(xmax = lo95_h2, xmin = up95_h2), size= 0.05) +
theme_cowplot(font_size= 8) +
scale_colour_manual(values= c('grey', colorBlindBlack8[c(8, 2, 4, 3)])) +
geom_vline(xintercept= 0, colour= colorBlindBlack8[1], linetype= 'dashed', size= 0.2, alpha= 0.6) +
geom_hline(yintercept= 0, colour= colorBlindBlack8[1], linetype= 'dashed', size= 0.2, alpha= 0.6) +
xlab('Effect size maternal \nnon-transmitted alleles, days') +
ylab('Effect size maternal genome, days') 
theme(legend.direction = "horizontal", legend.position = "bottom")
#scale_x_continuous(breaks = round(seq(-1.5, 3, by= 0.5), 1)) +
#  scale_y_continuous(breaks = round(seq(-1.5, 3, by= 0.5), 1))

ggsave(snakemake@output[[4]], plot= p1, width= 120, height= 60, units= 'mm', dpi= 300)
fwrite(d, snakemake@output[[5]], sep= '\t')
