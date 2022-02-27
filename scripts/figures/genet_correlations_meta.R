library(scales)
library("dplyr")
library("knitr")
library("tidyr")
library(cowplot)
library(ggrepel)
library("data.table")
library('showtext')
library(tidyverse)
library(fmsb)

colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

d= fread(snakemake@input[[1]])

d$p1= 'Gestational\nduration'
d$p2= with(d, ifelse(grepl('postTerm', p2), 'Post-term\ndelivery', ifelse(grepl('allPTD', p2), 'Preterm\ndelivery', 'GAnrm')))

d= filter(d, p2!= 'GAnrm')

p1= ggplot(d, aes(p2, rg, colour= p2)) +
  geom_point() +
geom_errorbar(aes(ymin= I(rg - 1.96*se) , ymax= (rg + 1.96 * se)), width=.2, position=position_dodge(.9)) +
theme_cowplot(font_size= 9) +
scale_fill_manual(values= colorBlindBlack8[c(8,3,2)], guide= 'none') +
scale_colour_manual(guide= 'none', values= colorBlindBlack8[c(8,3,2)]) +
xlab('Phenotype') +
ylab('Genetic correlation [95% CI]') +
theme(legend.position= 'none') +
ylim(pmin(-1, min(d$rg - 1.96*d$se)), pmax(1, max(d$rg + 1.96 * d$se))) +
geom_hline(yintercept= 0, linetype= 'dashed', colour= 'grey', size= 0.5)


ggsave(snakemake@output[[1]], plot= p1, width= 60, height= 80, units= 'mm', dpi= 300)

