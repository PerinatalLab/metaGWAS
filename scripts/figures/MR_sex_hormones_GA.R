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
x= fread(snakemake@input[[2]], select= c('ID', 'BETA', 'SE'))

mr= fread(snakemake@input[[3]])

d= inner_join(d,x, by= 'ID')
d= filter(d, !duplicated(ID))

d$BETA= with(d, ifelse(beta< 0, -1 * BETA, BETA))
d$beta= with(d, ifelse(beta< 0, -1 * beta, beta))


shbg= filter(d, trait== 'SHBG_fem_cluster')
testo= filter(d, trait== 'Testosterone_fem_cluster')

p1= ggplot(shbg, aes(beta, BETA), color= colorBlindBlack8[2]) +
geom_errorbarh(aes(xmin= beta - se, xmax= beta + se), size= 0.1, alpha= 0.7, color= colorBlindBlack8[2]) +
geom_errorbar(aes(ymin= BETA - SE, ymax= BETA + SE), size= 0.1, alpha= 0.7, color= colorBlindBlack8[2]) +
geom_point(size= 2, shape=21, stroke= 0.1, alpha= 0.7, fill= colorBlindBlack8[2]) +
xlab('Effect on SHBG (women), nmol/L') +
ylab('Effect on gestational duration, days') +
theme_cowplot(font_size= 8) +
geom_abline(intercept= 0, slope= filter(mr, method== 'IVW', trait== 'SHBG_fem_cluster') %>% pull(estimate), colour= '#d9d9d9') +
geom_abline(intercept= (filter(mr, method== '(intercept)', trait== 'SHBG_fem_cluster') %>% pull(estimate))[1], slope= filter(mr, method== 'MR-Egger', trait== 'SHBG_fem_cluster') %>% pull(estimate), colour= '#d9d9d9', linetype= 'dashed') +
geom_hline(yintercept= 0, size= 0.1) +
geom_vline(xintercept= 0, size= 0.1) +
theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks= element_blank(),
        panel.grid.major= element_line(colour= 'grey', size= 0.05))


p2= ggplot(testo, aes(beta, BETA), color= colorBlindBlack8[2]) +
geom_errorbarh(aes(xmin= beta - se, xmax= beta + se), size= 0.1, alpha= 0.7, color= colorBlindBlack8[2]) +
geom_errorbar(aes(ymin= BETA - SE, ymax= BETA + SE), size= 0.1, alpha= 0.7, color= colorBlindBlack8[2]) +
geom_point(size= 2, shape=21, stroke= 0.1, alpha= 0.7, fill= colorBlindBlack8[2]) +
xlab('Effect on testosterone (women), nmol/L') +
ylab('Effect on gestational duration, days') +
theme_cowplot(font_size= 8) +
geom_abline(intercept= 0, slope= filter(mr, method== 'IVW', trait== 'Testosterone_fem_cluster') %>% pull(estimate), colour= '#d9d9d9') +
geom_abline(intercept= (filter(mr, method== '(intercept)', trait== 'Testosterone_fem_cluster') %>% pull(estimate))[1], slope= filter(mr, method== 'MR-Egger', trait== 'Testosterone_fem_cluster') %>% pull(estimate), colour= '#d9d9d9', linetype= 'dashed') +
geom_hline(yintercept= 0, size= 0.1) +
geom_vline(xintercept= 0, size= 0.1) +
theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks= element_blank(),
        panel.grid.major= element_line(colour= 'grey', size= 0.05))

ggsave(snakemake@output[[1]], plot= p1, width= 70, height= 70, units= 'mm', dpi= 300)
ggsave(snakemake@output[[2]], plot= p2, width= 70, height= 70, units= 'mm', dpi= 300)
