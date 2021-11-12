library("dplyr")
library("knitr")
library("tidyr")
library(cowplot)
library(ggrepel)
library("data.table")
library('showtext')
options(warn=-1)


d= fread(snakemake@input[[1]], h= T)


colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

font_add("arial", "arial.ttf", bold= 'arial_bold.ttf')

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)



female_repr= c('breast', "cervix, uterine", 'endometrium', 'ovary', 'placenta', 'vagina', 'fallopian tube')
male_repr= c('ductus deferens', 'testis', 'seminal vesicle', 'prostate', 'epididymis')
muscle= c('smooth muscle', 'heart muscle', 'skeletal muscle')

d$organ= with(d, ifelse(tissue %in% female_repr, 'Female reproductive', ifelse(tissue %in% male_repr, 'Male reproductive', ifelse(tissue %in% muscle, 'Muscle', 'Others'))))



p1= ggplot(d, aes(-log10(MannW_pvalue), I(i_listmedian/ base_list_median), colour= organ)) +
geom_point(size= 1.5) +
theme_cowplot(font_size= 10) +
scale_colour_manual('Legend', values= c(colorBlindBlack8[c(3, 2, 8)], 'grey'), guide= 'none') +
geom_text_repel(data= filter(d, MannW_pvalue< 0.05), aes(label= tissue), fontface = 'bold') +
geom_vline(xintercept= -log10(0.05), colour= colorBlindBlack8[8], linetype= 'dashed', size= 0.2, alpha= 0.6) +
geom_vline(xintercept= -log10(0.05/nrow(d)), colour= colorBlindBlack8[8], linetype= 'dashed', size= 0.2, alpha= 0.6) +
ylab('Enrichment') +
xlab('-log10(pvalue)')


ggsave(snakemake@output[[1]], plot= p1, width= 120, height= 90, units= 'mm', dpi= 300)


p1= ggplot(d, aes(i_listmedian / base_list_median, tissue, colour= organ, size= -log10(MannW_pvalue))) +
geom_point() + 
theme_cowplot(font_size= 10) +
scale_colour_manual('Legend', values= c(colorBlindBlack8[c(3, 2, 8)], 'grey'), guide= 'none') +
scale_size_continuous(guide= 'none') +
xlab('Enrichment') +
ylab('') +
scale_x_continuous(limits= c(0, 2)) +
geom_vline(xintercept= seq(0, 2, 0.5), colour= 'grey', alpha= 0.4)



