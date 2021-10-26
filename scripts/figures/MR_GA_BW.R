library(MendelianRandomization)
library(data.table)
library(dplyr)
library("knitr")
library("tidyr")
library(cowplot)
library(ggrepel)
library('showtext')


colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

font_add("arial", "arial.ttf", bold= 'arial_bold.ttf')

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)

d= fread(snakemake@input[[1]])


top= fread(snakemake@input[[2]])
ids= pull(top, ID)
ids= c('3:156697097:A:G', '5:158058432:G:T', ids)

d= filter(d, ID %in% ids)
bw_m= fread(snakemake@input[[3]], select= c('ID', 'BETA', 'SE'))
names(bw_m)= c('ID', 'beta', 'se')

bw_f= fread(snakemake@input[[4]], select= c('ID', 'BETA', 'SE'))

names(bw_f)= c('ID', 'beta', 'se')

bw_m= inner_join(d, bw_m, by= 'ID')
bw_m$beta= with(bw_m, (ifelse(BETA<0, -1 * beta, beta)))
bw_m$BETA= with(bw_m, (ifelse(BETA<0, -1 * BETA, BETA)))

bw_f= inner_join(d, bw_f, by= 'ID')
bw_f$beta= with(bw_f, (ifelse(BETA<0, -1 * beta, beta)))
bw_f$BETA= with(bw_f, (ifelse(BETA<0, -1 * BETA, BETA)))

inputMR_m= mr_input(bx= bw_m$BETA, bxse= bw_m$SE, by= bw_m$beta, byse= bw_m$se)
mr_m= mr_allmethods(inputMR_m)$Values
names(mr_m)= c('method', 'estimate', 'se', 'lo95', 'up95', 'pvalue')

inputMR_f= mr_input(bx= bw_f$BETA, bxse= bw_f$SE, by= bw_f$beta, byse= bw_f$se)
mr_f= mr_allmethods(inputMR_f)$Values
names(mr_f)= c('method', 'estimate', 'se', 'lo95', 'up95', 'pvalue')


p1= ggplot(bw_m, aes(BETA, beta)) +
geom_errorbar(aes(ymin= beta - se, ymax= beta + se), alpha= 0.3, size= 0.3) +
geom_errorbarh(aes(xmin= BETA - SE, xmax= BETA + SE), alpha= 0.3, size= 0.3) +
geom_point(color= 'darkgrey', size= 0.5) +
xlab('Maternal effect \non gestational duration, days') +
ylab('Maternal only effect \non birth weight, z-score') +
theme_cowplot(font_size= 8) +
geom_abline(intercept= 0, slope= filter(mr_m, method== 'IVW') %>% pull(estimate), colour= colorBlindBlack8[2]) +
geom_abline(intercept= (filter(mr_m, method== '(intercept)') %>% pull(estimate))[1], slope= filter(mr_m, method== 'MR-Egger') %>% pull(estimate), colour= colorBlindBlack8[4]) +
geom_hline(yintercept= 0, size= 0.1) +
geom_vline(xintercept= 0, size= 0.1) +
theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks= element_blank(),
	panel.grid.major= element_line(colour= 'grey', size= 0.1))

ggsave(snakemake@output[[1]], plot= p1, width= 60, height= 60, units= 'mm', dpi= 300)

p1= ggplot(bw_f, aes(BETA, beta)) +
geom_errorbar(aes(ymin= beta - se, ymax= beta + se), alpha= 0.3, size= 0.3) +
geom_errorbarh(aes(xmin= BETA - SE, xmax= BETA + SE), alpha= 0.3, size= 0.3) +
geom_point(color= 'darkgrey', size= 0.5) +
xlab('Maternal effect \non gestational duration, days') +
ylab('Fetal only effect \non birth weight, z-score') +
theme_cowplot(font_size= 8) +
geom_abline(intercept= 0, slope= filter(mr_f, method== 'IVW') %>% pull(estimate), colour= colorBlindBlack8[2]) +
geom_abline(intercept= (filter(mr_f, method== '(intercept)') %>% pull(estimate))[1], slope= filter(mr_f, method== 'MR-Egger') %>% pull(estimate), colour= colorBlindBlack8[4]) +
geom_hline(yintercept= 0, size= 0.1) +
geom_vline(xintercept= 0, size= 0.1) +
theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks= element_blank(),
        panel.grid.major= element_line(colour= 'grey', size= 0.1))

ggsave(snakemake@output[[2]], plot= p1, width= 60, height= 60, units= 'mm', dpi= 300)


bw_f$trait= 'Fetal_BW'
bw_m$trait= 'Maternal_BW'

d= rbind(bw_f, bw_m)

mr_f$trait= 'Fetal_BW'
mr_m$trait= 'Maternal_BW'

fwrite(d, snakemake@output[[3]], sep= '\t')


d= rbind(mr_f, mr_m)
fwrite(d, snakemake@output[[4]], sep= '\t')

