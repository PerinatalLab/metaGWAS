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


top_ga= fread(snakemake@input[[1]])

top_ga= c(pull(top_ga, ID), '5:158058432:G:T', '3:156697097:A:G')


top_ptd= fread(snakemake@input[[2]])

top_ptd= pull(top_ptd, ID)

top= c(top_ga, top_ptd)

top= unique(top)

ga= fread(snakemake@input[[3]], select= c('ID', 'BETA', 'SE'))

ga= filter(ga, ID %in% top)

ptd= fread(snakemake@input[[4]], select= c('ID', 'BETA', 'SE'))

ptd= filter(ptd, ID %in% top) %>% select(ID, BETA, SE)

names(ptd)= c('ID', 'BETA_ptd', 'SE_ptd')

d= inner_join(ga, ptd, by= 'ID')

d$GWAS= with(d, ifelse(ID== '5:157895049:C:T', 'Both phenotypes', ifelse(ID %in% top_ptd, 'Preterm delivery', 'Gestational duration')))

p1= ggplot(d, aes(BETA, BETA_ptd, colour= GWAS, fill= GWAS)) +
geom_errorbarh(aes(xmin= BETA - SE, xmax= BETA + SE, colour= GWAS, fill= GWAS), size= 0.1, alpha= 0.7) +
geom_errorbar(aes(ymin= BETA_ptd - SE_ptd, ymax= BETA_ptd + SE_ptd, colour= GWAS, fill= GWAS),size= 0.1, alpha= 0.7) +
geom_point(size= 2, shape=21, stroke= 0.1, alpha= 0.7) +
scale_colour_manual(values= colorBlindBlack8[c(4, 2, 1)], guide= 'none') +
scale_fill_manual(values= colorBlindBlack8[c(4, 2, 1)], guide= 'none') +
xlab('Maternal effect on gestational duration, days') +
ylab('Maternal effect on preterm delivery, log(OR)') +
theme_cowplot(font_size= 8) +
geom_hline(yintercept= 0, size= 0.1) +
geom_vline(xintercept= 0, size= 0.1) +
theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks= element_blank(),
        panel.grid.major= element_line(colour= 'grey', size= 0.05))


ggsave(snakemake@output[[1]], plot= p1, width= 70, height= 70, units= 'mm', dpi= 300)

p1= ggplot(d, aes(BETA, BETA_ptd, colour= GWAS, fill= GWAS)) +
geom_errorbarh(aes(xmin= BETA - SE, xmax= BETA + SE, colour= GWAS, fill= GWAS), size= 0.1, alpha= 0.7) +
geom_errorbar(aes(ymin= BETA_ptd - SE_ptd, ymax= BETA_ptd + SE_ptd, colour= GWAS, fill= GWAS),size= 0.1, alpha= 0.7) +
geom_point(size= 2, shape=21, stroke= 0.1, alpha= 0.7) +
scale_colour_manual(values= colorBlindBlack8[c(4, 2, 1)], guide= 'none') +
scale_fill_manual(values= colorBlindBlack8[c(4, 2, 1)]) +
xlab('Maternal effect on gestational duration, days') +
ylab('Maternal effect on preterm delivery, log(OR)') +
theme_cowplot(font_size= 8) +
geom_hline(yintercept= 0, size= 0.1) +
geom_vline(xintercept= 0, size= 0.1) +
theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks= element_blank(),
        panel.grid.major= element_line(colour= 'grey', size= 0.05))

ggsave(snakemake@output[[2]], plot= p1, width= 70, height= 70, units= 'mm', dpi= 300)
