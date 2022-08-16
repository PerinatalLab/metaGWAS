library("dplyr")
library("knitr")
library("tidyr")
library(cowplot)
library(ggrepel)
library("data.table")
library('showtext')
options(warn=-1)

colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

font_add("arial", "arial.ttf", bold= 'arial_bold.ttf')

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)


d= fread(snakemake@input[[1]], head= T)

snp_list= c('rs2963463', 'rs201226733', 'rs56318008', 'rs4383453', 'rs200879388', 'rs2955117')
CHR_list= c(157895049, 115164770, 22470407, 123068359, 131300571, 127881613)
POS_list= c(5, 23, 1, 3, 23, 3)

df= data.frame(snp_list, CHR_list, POS_list)

d= inner_join(d, df, on= c('CHR', 'POS'))

x= fread(snakemake@input[[2]], head= T)

x= inner_join(x, df, on= c('CHR', 'POS'))

d= inner_join(d, x, by= c('CHR', 'POS'))

p1= ggplot(d, aes(BETA.y, BETA.x)) +
geom_errorbar(aes(ymin= BETA.x - SE.x, ymax= BETA.x + SE.x), alpha= 0.3, size= 0.3) +
geom_errorbarh(aes(xmin= BETA.y - SE.y, xmax= BETA.y + SE.y), alpha= 0.3, size= 0.3) +
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

fwrite()

