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


d= fread(snakemake@input[[1]])

d$term= with(d, ifelse(term== 'fetal_effect_PGS', 'Fetal', 'Maternal'))
d$outcome= gsub(' PGS', '', d$outcome)


p1= ggplot(d, aes(term, estimate, colour= term)) + 
geom_pointrange(aes(ymin= lo95, ymax= up95)) + 
facet_wrap(vars(outcome)) + 
scale_colour_manual(guide= 'none', values= colorBlindBlack8[c(2, 4)]) +
theme_cowplot(10) + 
geom_hline(yintercept= 0, colour= 'grey', size= 0.5, linetype= 'dashed') + 
theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) + 
ylab('Effect on gestational duration \ngenetic score (95% CI), days') +
xlab('Birth weight genetic score')


ggsave(snakemake@output[[1]], plot= p1, width= 180, height= 100, units= 'mm', dpi= 300)

